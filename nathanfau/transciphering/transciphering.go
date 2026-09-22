// Package transciphering wires the Algo1 IntRootBoot refresh into a homomorphic AES
// middle-round pipeline so an encrypted AES state can be advanced round by round under CKKS.
//
// One "middle" round (neither the first AddRoundKey-only round nor the last MixColumns-less
// round) is:
//
//	SubBytes -> Refresh(Algo1, ShiftRows at the pause) -> MixColumns -> AddRoundKey -> Cleaning
//
// Context.Round runs that sequence; each operation is exported on its own as well.
//
// The AES circuit runs in the conjugate-invariant (CI, real) context of the CtxSwitcher; the
// refresh round-trips each bit CI -> Std, packs k bits per ciphertext, bootstraps with Algo1,
// extracts the bits back, and returns to CI. The SubBytes version is selectable.
package transciphering

import (
	"fmt"
	"math"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
	"github.com/tuneinsight/lattigo/v6/nathanfau/algo1"
	"github.com/tuneinsight/lattigo/v6/nathanfau/bitbatching"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockpack"
	"github.com/tuneinsight/lattigo/v6/nathanfau/cleaning"
	"github.com/tuneinsight/lattigo/v6/nathanfau/convctx"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/nathanfau/params2"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

type Context struct {
	Params ckks.Parameters // Std bootstrapping / residual parameters
	Eval   *bootstrapping.Evaluator
	Sw     *convctx.CtxSwitcher
	AE     *aes.Evaluator // AES round functions, CI context (sw.EvalCI)

	EcdCI *ckks.Encoder
	EncCI *rlwe.Encryptor
	DecCI *rlwe.Decryptor

	K      int // bits packed per ciphertext (t = 2^K)
	Target int // Algo1 input level (S2C LevelQ)

	// Canon is the scale of the state in the AES circuit: where Refresh leaves it, and where the
	// state and the round keys are encrypted. It is the zone's scale when there is one, the default
	// scale otherwise.
	Canon rlwe.Scale
	// LastRefresh is where the last Refresh spent its time, stage by stage. The CSV writes it on the
	// Refresh row; nothing else reads it.
	LastRefresh RefreshTimings

	// Zone is the scale of the zone (params2.Shape.LogZone), zero when the chain has none. With a
	// zone, Refresh hops out of it onto algo1.EntryScale and Algo1 lands back on it.
	Zone rlwe.Scale

	Cfg Config // the variants this context runs on

	// Quiet silences the progress lines the pipeline prints, so a timed run measures the circuit
	// alone. The zero value keeps them.
	Quiet bool

	// SBoxExact runs SubBytes on the exactly aligned S-box (aes.SubByteExact). Temporary, to run a
	// full AES on it. It also lands the round cleaning on one fixed scale (see cleanFunc): the exact
	// S-box needs its 8 input bits on the same scale, which MixColumns does not leave them on.
	SBoxExact bool

	// CleanFixedScale lands the round cleaning on Canon rather than on its input's scale (see
	// cleanFunc). SBoxExact implies it.
	CleanFixedScale bool

	// RefreshCanon lands the refresh's bit extraction on Canon itself, so the bits come back on the
	// canonical scale exactly. Without it the extraction keeps the scale Algo1 hands it -- which
	// its squarings have pushed a few 2^-16 away from Canon -- and the refresh then RELABELS the
	// bits as Canon, which multiplies every 1 by that ratio.
	RefreshCanon bool

	// Levels that move with the cleaning depth, i.e. with the length of the chain, and with the
	// SlotsToCoeffs block: every prime it spends beyond the first lifts the pipeline by a level.
	SubBytesLv int // state level entering a round
	InitLv     int // rk0 and the encoded blocks: SubBytesLv + 1, one level for the XOR
	RefreshLv  int // level the refresh hands the state back at
	ARKLv      int // level AddRoundKey (or the fused XorClean) starts at

	// Levels the round keys have to be encrypted at, which follow Cfg.Place.
	ARKKeyLv  int // middle rounds
	LastKeyLv int // last round, no MixColumns

	// What the key material cost to build and how much of it there is, measured at construction.
	// The auxiliary modulus P moves both, so a P sweep reads them here.
	BtpKeyGen, SwKeyGen     time.Duration
	EvalSetup               time.Duration // the DFT matrices, which P does not move
	BtpKeyBytes, SwKeyBytes int
}

// Config selects the variants a Context runs on. The zero value is the pipeline default: the
// arithmetic XOR, the 2-level Cleaning, and AddRoundKey then Cleaning as two separate steps.
//
// Clean has to be fixed here rather than flipped later: its depth decides how many primes the
// chain carries, so it is baked into the parameters.
type Config struct {
	Xor   aes.XorKind
	Clean cleaning.Kind
	Place cleaning.Placement

	// CleanExtract runs the refresh on bitbatching.BitExtractClean instead of BitExtract: the
	// interpolation and the cleaning fused into one bivariate polynomial, so the bits come back
	// cleaned. One prime more on the chain, and an error quadratic in Algo1's output instead of
	// linear. Like Clean, it is baked into the parameters and cannot be flipped later.
	CleanExtract bool
}

// Extract is the extraction Cfg selects. Both have the same signature, the same output order and
// the same normalised scale and level, so the refresh does not care which one it holds.
func (c Config) Extract() func(ckks.Parameters, *ckks.Evaluator, *rlwe.Ciphertext, int) ([]*rlwe.Ciphertext, error) {
	if c.CleanExtract {
		return bitbatching.BitExtractClean
	}
	return bitbatching.BitExtract
}

// ExtractTo is Extract with the bits landed on a chosen scale.
func (c Config) ExtractTo() func(ckks.Parameters, *ckks.Evaluator, *rlwe.Ciphertext, int, rlwe.Scale) ([]*rlwe.Ciphertext, error) {
	if c.CleanExtract {
		return bitbatching.BitExtractCleanTo
	}
	return bitbatching.BitExtractTo
}

// ExtractLevels is how many primes that extraction spends, which the chain has to carry above the
// refresh.
func (c Config) ExtractLevels(k int) int {
	if c.CleanExtract {
		return k + 1
	}
	return k
}

// The levels below are those of the DEFAULT chain, whose SlotsToCoeffs spends a single prime; a
// context carries its own in SubBytesLv, InitLv, RefreshLv and ARKLv, which is what the pipeline
// reads. These stay as the reference the tests are written against.
const (
	SubBytesLevel = 5 // state level entering a round; below the cleaning block, so it never moves
	InitLevel     = 6 // rk0 and the encoded blocks: SubBytesLevel + 1, one level for the XOR

	// DefaultCleanDepth is cleaning.Basic's, i.e. the chain Config's zero value builds.
	DefaultCleanDepth = 2
	RefreshLevel      = SubBytesLevel + 4 + DefaultCleanDepth // 11
	ARKLevel          = RefreshLevel - 3                      // 8
)

// RefreshLevelFor and ARKLevelFor are RefreshLevel and ARKLevel on a chain built for a cleaning
// polynomial of the given depth. Between the refresh and the next round the state pays MixColumns
// (3), then the XOR (1) and the cleaning (cleanDepth), whatever the placement.
func RefreshLevelFor(cleanDepth int) int { return SubBytesLevel + 4 + cleanDepth }
func ARKLevelFor(cleanDepth int) int     { return RefreshLevelFor(cleanDepth) - 3 }

// NewContext builds the whole machinery on the algo1 moduli chain, in its default configuration.
// The AES state lives in the CI ring (sw.CiP), the bootstrap runs in the Std ring.
func NewContext(logN, k int) (*Context, error) {
	return NewContextWith(logN, k, Config{})
}

// NewContextWith is NewContext on the variants cfg selects. The chain is built for cfg.Clean's
// depth, so RefreshLv and ARKLv follow it.
func NewContextWith(logN, k int, cfg Config) (*Context, error) {
	return NewContextWith2(logN, k, cfg, params2.Shape{})
}

// NewContextWith2 is NewContextWith on a chosen parameter shape (see params2.Shape). P leaves Q
// untouched and moves only the key-switching, where the SlotsToCoeffs block and q0 reshape the
// chain itself, so the context's own levels follow them. A zone keeps the levels and moves the
// scale instead: see Context.Zone.
func NewContextWith2(logN, k int, cfg Config, sh params2.Shape) (*Context, error) {
	depth := cfg.Clean.Depth()
	params, btpParams, err := params2.TranscipheringParamsWith(logN, k, depth, cfg.ExtractLevels(k), sh)
	if err != nil {
		return nil, fmt.Errorf("TranscipheringParams: %w", err)
	}

	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()

	fmt.Println("BTS KeyGen ...")
	t0 := time.Now()
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		return nil, fmt.Errorf("GenEvaluationKeys: %w", err)
	}
	btpKeyGen, btpKeyBytes := time.Since(t0), evk.BinarySize()
	t1 := time.Now()
	eval, err := bootstrapping.NewEvaluator(btpParams, evk)
	if err != nil {
		return nil, fmt.Errorf("bootstrapping NewEvaluator: %w", err)
	}
	// NewEvaluator is not key generation: it encodes the two DFT matrices. It costs real time, so it
	// is timed apart.
	evalSetup := time.Since(t1)
	fmt.Printf("Done !(%s, dont %s de matrices DFT)\n", time.Since(t0).Round(time.Millisecond), evalSetup.Round(time.Millisecond))

	fmt.Println("CtxSwitcher KeyGen ...")
	t0 = time.Now()
	sw, err := convctx.NewCtxSwitcher(params, sk)
	if err != nil {
		return nil, fmt.Errorf("NewCtxSwitcher: %w", err)
	}
	swKeyGen := time.Since(t0)
	fmt.Printf("Done !(%s)\n", swKeyGen.Round(time.Millisecond))

	ae := aes.NewEvaluator(sw.EvalCI)

	c := &Context{
		Params: params,
		Eval:   eval,
		Sw:     sw,
		AE:     ae,
		EcdCI:  ckks.NewEncoder(sw.CiP),
		EncCI:  rlwe.NewEncryptor(sw.CiP, sw.SkCI),
		DecCI:  rlwe.NewDecryptor(sw.CiP, sw.SkCI),
		K:      k,
		Target: eval.SlotsToCoeffsParameters.LevelQ,
		Canon:  sw.CiP.DefaultScale(),

		Cfg: cfg,

		BtpKeyGen:   btpKeyGen,
		SwKeyGen:    swKeyGen,
		EvalSetup:   evalSetup,
		BtpKeyBytes: btpKeyBytes,
		SwKeyBytes:  sw.KeyBytes,
	}

	// Target is the SlotsToCoeffs prime count: one is the default chain, each extra one lifts
	// everything above it by a level.
	lift := c.Target - 1
	c.SubBytesLv = SubBytesLevel + lift
	c.InitLv = InitLevel + lift
	c.RefreshLv = RefreshLevelFor(depth) + lift
	c.ARKLv = ARKLevelFor(depth) + lift
	// CleanOne is the only placement that sends the key straight to the XOR, i.e. below h(state).
	c.ARKKeyLv = c.ARKLv - cfg.Place.YDrop(depth)
	c.LastKeyLv = c.RefreshLv - cfg.Place.YDrop(depth)

	// With a zone the state lives at its scale, so Canon follows it. Without one nothing hops, and
	// the pipeline runs exactly as it did before zones existed.
	if sh.HasZone() {
		c.Zone = sh.ZoneScale()
		c.Canon = c.Zone
	}

	// Debug decoding contexts: Std (residual) and CI (AES circuit).
	debug.EncStd, debug.DecStd, debug.ParamsStd = ckks.NewEncoder(params), rlwe.NewDecryptor(params, sk), params
	debug.EncCI, debug.DecCI, debug.ParamsCI = c.EcdCI, c.DecCI, sw.CiP

	return c, nil
}

// MemItem is one object a run keeps alive, and what it weighs.
type MemItem struct {
	Name  string
	Count int
	Bytes int
}

// MemBreakdown lists the big objects the context holds, weighed on the objects themselves (their
// serialized size, which for lattigo's polynomials is their memory to a few bytes). What it does not
// list -- evaluator buffers, encoders, ring tables, Go overhead -- is the gap to the live heap.
func (c *Context) MemBreakdown() []MemItem {
	var items []MemItem
	add := func(name string, count, bytes int) {
		if count > 0 {
			items = append(items, MemItem{name, count, bytes})
		}
	}
	evk := c.Eval.EvaluationKeys
	n, b := evkSize(evk.EvkN1ToN2, evk.EvkN2ToN1, evk.EvkRealToCmplx, evk.EvkCmplxToReal)
	add("btp ring-switch keys", n, b)
	n, b = evkSize(evk.EvkDenseToSparse, evk.EvkSparseToDense)
	add("btp sparse-secret keys", n, b)
	if ks := evk.MemEvaluationKeySet; ks != nil {
		if ks.RelinearizationKey != nil {
			add("btp relinearization key", 1, ks.RelinearizationKey.BinarySize())
		}
		b = 0
		for _, gk := range ks.GaloisKeys {
			b += gk.BinarySize()
		}
		add("btp Galois keys", len(ks.GaloisKeys), b)
	}
	n, b = dftSize(c.Eval.S2CDFTMatrix)
	add("DFT SlotsToCoeffs (QP)", n, b)
	n, b = dftSize(c.Eval.C2SDFTMatrix)
	add("DFT CoeffsToSlots (QP)", n, b)
	add("switcher keys (deg 2N)", 3, c.Sw.KeyBytes)
	n, b = c.Sw.MaskBytes()
	add("switcher masks (deg 2N)", n, b)
	return items
}

// DFTSize is the two DFT matrices: plaintext diagonals and bytes, SlotsToCoeffs then CoeffsToSlots.
func (c *Context) DFTSize() (s2cPts, s2cBytes, c2sPts, c2sBytes int) {
	s2cPts, s2cBytes = dftSize(c.Eval.S2CDFTMatrix)
	c2sPts, c2sBytes = dftSize(c.Eval.C2SDFTMatrix)
	return
}

func evkSize(keys ...*rlwe.EvaluationKey) (n, bytes int) {
	for _, k := range keys {
		if k != nil {
			n++
			bytes += k.BinarySize()
		}
	}
	return n, bytes
}

// dftSize counts a DFT matrix's plaintext diagonals and their size.
func dftSize(m dft.Matrix) (plaintexts, bytes int) {
	for _, lt := range m.Matrices {
		for _, p := range lt.Vec {
			plaintexts++
			bytes += p.BinarySize()
		}
	}
	return plaintexts, bytes
}

// logf prints a progress line unless the context is Quiet.
func (c *Context) logf(format string, a ...any) {
	if !c.Quiet {
		fmt.Printf(format, a...)
	}
}

// FirstRound is AES round 0: it XORs rk0 into the AES blocks. The blocks arrive in the CLEAR,
// so the 'server' encodes them itself and the XOR is ciphertext against plaintext.
// rk0 must be at InitLevel.
func (c *Context) FirstRound(blocks [][16]byte, rk0 blockpack.Packed) (blockpack.Packed, error) {
	if l := rk0[0][0].Level(); l != c.InitLv {
		return blockpack.Packed{}, fmt.Errorf("FirstRound: rk0 at level %d, want InitLevel %d", l, c.InitLv)
	}
	xor := c.AE.PlainOf(c.Cfg.Xor)

	var out blockpack.Packed
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			x, err := xor(rk0[g][b], blockpack.SlotVec(c.Sw.CiP, blocks, g, b))
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("FirstRound xor [%d][%d]: %w", g, b, err)
			}
			out[g][b] = x
		}
	}
	return out, nil
}

// SubBytes applies the selected SubByte version to the whole state, i.e. to all 64 ciphertexts.
func (c *Context) SubBytes(st blockpack.Packed, version int) (blockpack.Packed, error) {
	sbox := c.AE.SubByte
	if c.SBoxExact {
		sbox = c.AE.SubByteExact
	}
	var out blockpack.Packed
	for g := 0; g < 8; g++ {
		t0 := time.Now()
		ob, err := sbox(st[g], version)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("SubBytes group %d (v%d): %w", g, version, err)
		}
		c.logf("[SubBytes v%d] group %2d/8 done (%s)\n", version, g+1, time.Since(t0).Round(time.Millisecond))
		out[g] = ob
	}
	return out, nil
}

// RefreshTimings is where the last Refresh spent its time. It embeds algo1's own breakdown and adds
// what the refresh does around it -- ShiftRows, the bit extraction, and the conversions that are
// not Algorithm 1's. Reset at the top of every Refresh, so it always describes the LAST one.
//
// The buckets are meant to ACCOUNT FOR the whole refresh: Other absorbs BitPack, ScaleDown, ModUp
// and the recombinations, so STC + CTS + EvalMod + Extract + Conv + SR + Other comes back to the
// refresh's own duration, give or take the scheduler. A breakdown that does not add up hides
// exactly what one wants to see.
type RefreshTimings struct {
	algo1.Timings
	BitPack   time.Duration // packing the 4-bit nibbles in CI, before the conversion out
	SR        time.Duration // ShiftRows at the Algo1 pause
	Recombine time.Duration // CombineReIm on the 64 extracted bits, after the pause
}

// zoneScaleTol is how close, in bits, Refresh's output must land on the zone's scale before it is
// labelled with it: off by more, the label would change the value.
const zoneScaleTol = 30

// inZone reports whether the context runs a zone (see Context.Zone).
func (c *Context) inZone() bool { return c.Zone.Value.Sign() != 0 }

// Refresh refreshes the whole state and applies ShiftRows at the Algo1 pause. With a zone, the
// conversion to Std is also the hop out of it, onto the scale Algo1 wants, and Algo1 lands back on
// the zone's scale (algo1.ExtractTo).
func (c *Context) Refresh(st blockpack.Packed) (blockpack.Packed, error) {
	eval := c.Eval
	ck := eval.Evaluator
	entry := algo1.EntryScale(eval)
	c.LastRefresh = RefreshTimings{} // ce chrono decrit LE refresh en cours, pas la somme du run
	tm := &c.LastRefresh

	// 1-2. BitPack the 4-bit nibbles in CI, then convert to Std
	t0 := time.Now()
	var packed [16]*rlwe.Ciphertext
	for g := 0; g < 8; g++ {
		tPack := time.Now()
		lo, err := bitbatching.BitPack(c.Sw.EvalCI, st[g][0:4])
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh BitPack low g=%d: %w", g, err)
		}
		hi, err := bitbatching.BitPack(c.Sw.EvalCI, st[g][4:8])
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh BitPack high g=%d: %w", g, err)
		}
		tm.BitPack += time.Since(tPack)
		for i, ciNib := range [2]*rlwe.Ciphertext{lo, hi} {
			// Algo1 veut son entree a EntryScale, zone ou pas. ScaleDown ne corrige l'ecart qu'a un
			// facteur ENTIER pres ; ce qui reste atteint le mod-1 comme une erreur RELATIVE sur le
			// message, qu'Algo1, contrairement a EvalMod, ne sait pas rediviser ensuite. Cette
			// conversion est gratuite (elle multiplie par une constante), donc la sauter ne faisait
			// economiser rien et coutait 6,6 bits sur la sortie d'Algo1 : mesure a logN 16, sans
			// zone, l'erreur d'Algo1 valait 2^-5,9 contre 2^-12,5 avec une zone, ce qui tuait
			// l'extraction `half` (lineaire en cette erreur) des le tour 3 et laissait passer
			// `clean` (quadratique) a exactement le double de bits.
			tConv := time.Now()
			s, err := c.Sw.CIToStandardTo(ciNib, entry)
			tm.Conv += time.Since(tConv)
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("refresh CIToStandard g=%d nibble %d: %w", g, i, err)
			}
			packed[2*g+i] = s
		}
	}
	c.logf("[Algo1] BitPack + CI->Std, 16 packets done (%s)\n", time.Since(t0).Round(time.Millisecond))

	// 3. algo1.Extract per packet (drop to the Algo1 input level first).
	var reals, imags [16]*rlwe.Ciphertext
	for p := 0; p < 16; p++ {
		tPkt := time.Now()
		if d := packed[p].Level() - c.Target; d > 0 {
			ck.DropLevel(packed[p], d)
		}
		rr, ii, err := algo1.ExtractToTimed(eval, c.Sw, packed[p], c.K, c.Zone, &tm.Timings)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh Extract packet %d: %w", p, err)
		}
		reals[p], imags[p] = rr, ii
		c.logf("[Algo1] packet %2d/16: STC + ModUp + CTS + cos/sin + Re/Im split done (%s)\n", p+1, time.Since(tPkt).Round(time.Millisecond))
	}

	// 4. ShiftRows at the pause: pure pointer moves on the 32 packed nibbles.
	tSR := time.Now()
	re, im := c.AE.ShiftRows(reals, imags)
	tm.SR += time.Since(tSR)

	// 5. algo1.Resume: the double-angle squarings, on all 32 nibbles (square is pointwise, commutes
	//    with the ShiftRows permutation).
	t0 = time.Now()
	for k := 0; k < 16; k++ {
		if err := algo1.ResumeTimed(eval, re[k], im[k], &tm.Timings); err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh Resume %d: %w", k, err)
		}
	}
	c.logf("[Algo1] ShiftRows + squarings, 32 nibbles done (%s)\n", time.Since(t0).Round(time.Millisecond))

	// 6. Per output group, the bit extraction of its 4 nibbles (its two bytes), then CombineReIm (real
	//    byte + i*imag byte), Std -> CI and the canonical scale on its 8 bits.
	// The extraction lands on Canon under RefreshCanon, on its input's scale otherwise. The Re/Im
	// recombination (a multiplication by i) and StandardToCI both keep the scale, so what reaches
	// the relabel below is exactly what the extraction produced.
	extractTo := c.Cfg.ExtractTo()
	extract := func(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext, k int) ([]*rlwe.Ciphertext, error) {
		if c.RefreshCanon {
			return extractTo(params, eval, ct, k, c.Canon)
		}
		return extractTo(params, eval, ct, k, ct.Scale)
	}
	relabelGap := math.Inf(1) // the largest scale gap the relabel below overwrites, as -log2
	var out blockpack.Packed
	for j := 0; j < 8; j++ {
		tExt := time.Now()
		aLo, err := extract(c.Params, ck, re[2*j], c.K)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh extract re low j=%d: %w", j, err)
		}
		aHi, err := extract(c.Params, ck, re[2*j+1], c.K)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh extract re high j=%d: %w", j, err)
		}
		bLo, err := extract(c.Params, ck, im[2*j], c.K)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh extract im low j=%d: %w", j, err)
		}
		bHi, err := extract(c.Params, ck, im[2*j+1], c.K)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh extract im high j=%d: %w", j, err)
		}
		dExt := time.Since(tExt)
		tm.Extract += dExt
		aBits := append(append([]*rlwe.Ciphertext{}, aLo...), aHi...) // 8 bits of the real byte
		bBits := append(append([]*rlwe.Ciphertext{}, bLo...), bHi...) // 8 bits of the imag byte
		tConv := time.Now()
		for beta := 0; beta < 8; beta++ {
			tRe := time.Now()
			z, err := utils.CombineReIm(ck, aBits[beta], bBits[beta])
			tm.Recombine += time.Since(tRe)
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("refresh CombineReIm j=%d beta=%d: %w", j, beta, err)
			}
			tStd := time.Now()
			ci, err := c.Sw.StandardToCI(z)
			tm.Conv += time.Since(tStd)
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("refresh StandardToCI j=%d beta=%d: %w", j, beta, err)
			}
			gap := ci.Scale.Log2Delta(c.Canon)
			if c.RefreshCanon && gap < 40 {
				return blockpack.Packed{}, fmt.Errorf("refresh j=%d beta=%d: extraction landed 2^-%.1f off Canon", j, beta, gap)
			}
			if c.inZone() && gap < zoneScaleTol {
				return blockpack.Packed{}, fmt.Errorf("refresh j=%d beta=%d: scale 2^%.6f is 2^-%.1f off the zone's 2^%.6f, want below 2^-%d",
					j, beta, ci.Scale.Log2(), gap, c.Canon.Log2(), zoneScaleTol)
			}
			relabelGap = math.Min(relabelGap, gap)
			ci.Scale = c.Canon
			out[j][beta] = ci
		}
		c.logf("[Algo1] group %d/8: bit extraction, 4 nibbles (%s) + Re/Im recombination and Std->CI, 8 bits (%s)\n",
			j+1, dExt.Round(time.Millisecond), time.Since(tConv).Round(time.Millisecond))
	}
	if math.IsInf(relabelGap, 1) {
		c.logf("[Algo1] bits on Canon exactly, nothing relabelled\n")
	} else {
		c.logf("[Algo1] relabel to Canon: true scale off by up to 2^-%.1f (a relative error of that size on every 1)\n", relabelGap)
	}
	return out, nil
}

// AddRoundKey XORs the round key into the state, bit by bit. rk must already be at ARKLevel
func (c *Context) AddRoundKey(st, rk blockpack.Packed) (blockpack.Packed, error) {
	if kl, sl := rk[0][0].Level(), st[0][0].Level(); kl != sl {
		return blockpack.Packed{}, fmt.Errorf("AddRoundKey: round key at level %d, state at level %d", kl, sl)
	}
	xor := c.AE.Of(c.Cfg.Xor)
	var out blockpack.Packed
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			x, err := xor(st[g][b], rk[g][b])
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("AddRoundKey xor [%d][%d]: %w", g, b, err)
			}
			out[g][b] = x
		}
	}
	return out, nil
}

// Clean pulls every bit of the state back onto 0 and 1 after the round has spread them, with the
// polynomial Cfg.Clean names (SmootherCleaning by default, see discussion in []).
func (c *Context) Clean(st blockpack.Packed) (blockpack.Packed, error) {
	clean := c.cleanFunc()
	var out blockpack.Packed
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			cc, err := clean(c.Sw.CiP, c.Sw.EvalCI, st[g][b])
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("Clean [%d][%d]: %w", g, b, err)
			}
			out[g][b] = cc
		}
	}
	return out, nil
}

// cleanFunc is the round cleaning. The MixColumns trees have 5 or 7 leaves depending on the bit, so
// the 8 bits of a byte reach the cleaning on scales a few 2^-17 apart, and a cleaning that keeps
// its input's scale hands them on to the next SubBytes. The exact S-box refuses that, so under
// SBoxExact (or CleanFixedScale) the cleaning lands every bit on Canon instead -- the default scale,
// or the zone's -- no level, the polynomial
// evaluator picks its constants for the target.
func (c *Context) cleanFunc() cleaning.CleanFunc {
	if c.CleanFixedScale || c.SBoxExact {
		return c.Cfg.Clean.FuncAt(c.Canon)
	}
	return c.Cfg.Clean.Func()
}

// XorClean is AddRoundKey and Clean in one step, for the same depth: it cleans the state and the
// round key (CleanBoth) or the state only (CleanOne) BEFORE XORing them, instead of XORing first
// and cleaning the result. Cleaning before the XOR keeps the polynomial inside its basin and does
// not let the XOR amplify what the cleaning just contracted. rk must be at ARKKeyLv.
func (c *Context) XorClean(st, rk blockpack.Packed) (blockpack.Packed, error) {
	drop := c.Cfg.Place.YDrop(c.Cfg.Clean.Depth())
	if kl, sl := rk[0][0].Level(), st[0][0].Level(); kl != sl-drop {
		return blockpack.Packed{}, fmt.Errorf("XorClean (%s): round key at level %d, state at level %d, want key at %d",
			c.Cfg.Place, kl, sl, sl-drop)
	}
	route, xor, clean := c.Cfg.Place.Func(), c.AE.Of(c.Cfg.Xor), c.cleanFunc()

	var out blockpack.Packed
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			x, err := route(c.Sw.CiP, c.Sw.EvalCI, clean, xor, st[g][b], rk[g][b])
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("XorClean [%d][%d]: %w", g, b, err)
			}
			out[g][b] = x
		}
	}
	return out, nil
}

// arkThenClean runs the tail of a round: CleanAfter as the two historical steps, so the caller
// still gets a trace between the XOR and the polynomial, the other two as one fused XorClean.
func (c *Context) arkThenClean(st, rk blockpack.Packed, after RoundStep) (blockpack.Packed, error) {
	if c.Cfg.Place != cleaning.CleanAfter {
		c.logf("---- XorClean ----\n")
		t0 := time.Now()
		out, err := c.XorClean(st, rk)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("XorClean: %w", err)
		}
		after("XorClean", time.Since(t0), out)
		return out, nil
	}

	c.logf("---- AddRoundKey ----\n")
	t0 := time.Now()
	st, err := c.AddRoundKey(st, rk)
	if err != nil {
		return blockpack.Packed{}, fmt.Errorf("AddRoundKey: %w", err)
	}
	after("AddRoundKey", time.Since(t0), st)

	c.logf("---- Cleaning ----\n")
	t0 = time.Now()
	if st, err = c.Clean(st); err != nil {
		return blockpack.Packed{}, fmt.Errorf("Cleaning: %w", err)
	}
	after("Cleaning", time.Since(t0), st)
	return st, nil
}

// RoundStep is called after each operation of a middle round: its name, its duration, the state it
// produced. That is where a caller hooks its timings and its oracle checks; nil skips it.
type RoundStep func(name string, dur time.Duration, st blockpack.Packed)

// Round advances the packed state by one AES middle round: SubBytes -> Refresh -> MixColumns ->
// AddRoundKey -> Cleaning. rk is the round key at ARKKeyLv, version the SubByte variant (1..3).
func (c *Context) Round(st, rk blockpack.Packed, version int, after RoundStep) (blockpack.Packed, error) {
	if after == nil {
		after = func(string, time.Duration, blockpack.Packed) {}
	}
	var err error

	c.logf("---- SubBytes ----\n")
	t0 := time.Now()
	if st, err = c.SubBytes(st, version); err != nil {
		return blockpack.Packed{}, fmt.Errorf("Round SubBytes: %w", err)
	}
	after("SubBytes", time.Since(t0), st)

	c.logf("---- Refresh ----\n")
	t0 = time.Now()
	if st, err = c.Refresh(st); err != nil {
		return blockpack.Packed{}, fmt.Errorf("Round Refresh: %w", err)
	}
	after("Refresh", time.Since(t0), st)

	c.logf("---- MixColumns ----\n")
	t0 = time.Now()
	if st, err = c.AE.MixColumnsWith(c.AE.Of(c.Cfg.Xor), st); err != nil {
		return blockpack.Packed{}, fmt.Errorf("Round MixColumns: %w", err)
	}
	after("MixColumns", time.Since(t0), st)

	if st, err = c.arkThenClean(st, rk, after); err != nil {
		return blockpack.Packed{}, fmt.Errorf("Round: %w", err)
	}
	return st, nil
}

// LastRoundV1 is AES round 10: Round without MixColumns. ShiftRows still has to happen, and it stays
// inside the Refresh -- on a packed state, swapping the Re/Im halves of a group would take a
// half-slot rotation. Skipping MixColumns means the state reaches
// AddRoundKey at RefreshLevel, so rk must be encrypted there and NOT at ARKLevel.
func (c *Context) LastRoundV1(st, rk blockpack.Packed, version int, after RoundStep) (blockpack.Packed, error) {
	if after == nil {
		after = func(string, time.Duration, blockpack.Packed) {}
	}
	var err error

	c.logf("---- SubBytes ----\n")
	t0 := time.Now()
	if st, err = c.SubBytes(st, version); err != nil {
		return blockpack.Packed{}, fmt.Errorf("LastRoundV1 SubBytes: %w", err)
	}
	after("SubBytes", time.Since(t0), st)

	c.logf("---- Refresh ----\n")
	t0 = time.Now()
	if st, err = c.Refresh(st); err != nil {
		return blockpack.Packed{}, fmt.Errorf("LastRoundV1 Refresh: %w", err)
	}
	after("Refresh", time.Since(t0), st)

	if st, err = c.arkThenClean(st, rk, after); err != nil {
		return blockpack.Packed{}, fmt.Errorf("LastRoundV1: %w", err)
	}
	return st, nil
}
