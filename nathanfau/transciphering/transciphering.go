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
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
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

	K      int        // bits packed per ciphertext (t = 2^K)
	Target int        // Algo1 input level (S2C LevelQ)
	Canon  rlwe.Scale // canonical output scale after refresh
}

const (
	SubBytesLevel = 5  // state level entering a round
	InitLevel     = 6  // rk0 and the encoded blocks: SubBytesLevel + 1, the square of the XOR
	RefreshLevel  = 12 // level the refresh hands the state back at
	ARKLevel      = 9  // RefreshLevel - 3, the depth of MixColumns' XOR tree
)

// NewContext builds the whole machinery on the algo1 moduli chain.
// The AES state lives in the CI ring (sw.CiP)
// The bootstrap runs in the Std ring
func NewContext(logN, k int) (*Context, error) {
	params, btpParams, err := params2.TranscipheringParams(logN, k)
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
	eval, err := bootstrapping.NewEvaluator(btpParams, evk)
	if err != nil {
		return nil, fmt.Errorf("bootstrapping NewEvaluator: %w", err)
	}
	fmt.Printf("Done !(%s)\n", time.Since(t0).Round(time.Millisecond))

	fmt.Println("CtxSwitcher KeyGen ...")
	t0 = time.Now()
	sw, err := convctx.NewCtxSwitcher(params, sk)
	if err != nil {
		return nil, fmt.Errorf("NewCtxSwitcher: %w", err)
	}
	fmt.Printf("Done !(%s)\n", time.Since(t0).Round(time.Millisecond))

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
	}

	// Debug decoding contexts: Std (residual) and CI (AES circuit).
	debug.EncStd, debug.DecStd, debug.ParamsStd = ckks.NewEncoder(params), rlwe.NewDecryptor(params, sk), params
	debug.EncCI, debug.DecCI, debug.ParamsCI = c.EcdCI, c.DecCI, sw.CiP

	return c, nil
}

// FirstRound is AES round 0: it XORs rk0 into the AES blocks. The blocks arrive in the CLEAR,
// so the 'server' encodes them itself and the XOR is ciphertext against plaintext.
// rk0 must be at InitLevel.
func (c *Context) FirstRound(blocks [][16]byte, rk0 blockpack.Packed) (blockpack.Packed, error) {
	if l := rk0[0][0].Level(); l != InitLevel {
		return blockpack.Packed{}, fmt.Errorf("FirstRound: rk0 at level %d, want InitLevel %d", l, InitLevel)
	}
	eval := c.Sw.EvalCI

	var out blockpack.Packed
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			pt := ckks.NewPlaintext(c.Sw.CiP, InitLevel)
			pt.Scale = rk0[g][b].Scale
			if err := c.EcdCI.Encode(blockpack.SlotVec(c.Sw.CiP, blocks, g, b), pt); err != nil {
				return blockpack.Packed{}, fmt.Errorf("FirstRound encode [%d][%d]: %w", g, b, err)
			}
			x, err := eval.SubNew(rk0[g][b], pt)
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("FirstRound sub [%d][%d]: %w", g, b, err)
			}
			if err := eval.MulRelin(x, x, x); err != nil {
				return blockpack.Packed{}, fmt.Errorf("FirstRound square [%d][%d]: %w", g, b, err)
			}
			if err := eval.Rescale(x, x); err != nil {
				return blockpack.Packed{}, fmt.Errorf("FirstRound rescale [%d][%d]: %w", g, b, err)
			}
			out[g][b] = x
		}
	}
	return out, nil
}

// SubBytes applies the selected SubByte version to the whole state, i.e. to all 64 ciphertexts.
func (c *Context) SubBytes(st blockpack.Packed, version int) (blockpack.Packed, error) {
	var out blockpack.Packed
	for g := 0; g < 8; g++ {
		t0 := time.Now()
		ob, err := c.AE.SubByte(st[g], version)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("SubBytes group %d (v%d): %w", g, version, err)
		}
		fmt.Printf("[SubBytes v%d] group %2d/8 done (%s)\n", version, g+1, time.Since(t0).Round(time.Millisecond))
		out[g] = ob
	}
	return out, nil
}

// Refresh refreshes the whole state and applies ShiftRows at the Algo1 pause
func (c *Context) Refresh(st blockpack.Packed) (blockpack.Packed, error) {
	eval := c.Eval
	ck := eval.Evaluator

	// 1. CI -> Std, per bit.
	var stdFold [8][8]*rlwe.Ciphertext
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			s, err := c.Sw.CIToStandard(st[g][b])
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("refresh CIToStandard [%d][%d]: %w", g, b, err)
			}
			stdFold[g][b] = s
		}
	}

	// 2. BitPack 4-bit nibbles: packet 2g = low nibble (bits 0..3), 2g+1 = high nibble (bits 4..7);
	//    each packet is byte_g(nibble) + i*byte_{g+8}(nibble).
	var packed [16]*rlwe.Ciphertext
	for g := 0; g < 8; g++ {
		lo, err := bitbatching.BitPack(ck, stdFold[g][0:4])
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh BitPack low g=%d: %w", g, err)
		}
		hi, err := bitbatching.BitPack(ck, stdFold[g][4:8])
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh BitPack high g=%d: %w", g, err)
		}
		packed[2*g], packed[2*g+1] = lo, hi
	}

	// 3. algo1.Extract per packet (drop to the Algo1 input level first).
	var reals, imags [16]*rlwe.Ciphertext
	for p := 0; p < 16; p++ {
		tPkt := time.Now()
		if d := packed[p].Level() - c.Target; d > 0 {
			ck.DropLevel(packed[p], d)
		}
		rr, ii, err := algo1.Extract(eval, c.Sw, packed[p], c.K)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh Extract packet %d: %w", p, err)
		}
		reals[p], imags[p] = rr, ii
		fmt.Printf("[refresh] extract packet %2d/16 done (%s)\n", p+1, time.Since(tPkt).Round(time.Millisecond))
	}

	// 4. ShiftRows at the pause: pure pointer moves on the 32 packed nibbles.
	re, im := c.AE.ShiftRows(reals, imags)

	// 5. algo1.Resume: the double-angle squarings, on all 32 nibbles (square is pointwise, commutes
	//    with the ShiftRows permutation).
	for k := 0; k < 16; k++ {
		if err := algo1.Resume(eval, re[k], im[k]); err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh Resume %d: %w", k, err)
		}
	}

	// 6. BitExtract each nibble, CombineReIm (real byte + i*imag byte), Std -> CI, canonical scale.
	var out blockpack.Packed
	for j := 0; j < 8; j++ {
		aLo, err := bitbatching.BitExtract(c.Params, ck, re[2*j], c.K)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh BitExtract re low j=%d: %w", j, err)
		}
		aHi, err := bitbatching.BitExtract(c.Params, ck, re[2*j+1], c.K)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh BitExtract re high j=%d: %w", j, err)
		}
		bLo, err := bitbatching.BitExtract(c.Params, ck, im[2*j], c.K)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh BitExtract im low j=%d: %w", j, err)
		}
		bHi, err := bitbatching.BitExtract(c.Params, ck, im[2*j+1], c.K)
		if err != nil {
			return blockpack.Packed{}, fmt.Errorf("refresh BitExtract im high j=%d: %w", j, err)
		}
		aBits := append(append([]*rlwe.Ciphertext{}, aLo...), aHi...) // 8 bits of the real byte
		bBits := append(append([]*rlwe.Ciphertext{}, bLo...), bHi...) // 8 bits of the imag byte
		for beta := 0; beta < 8; beta++ {
			z, err := utils.CombineReIm(ck, aBits[beta], bBits[beta])
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("refresh CombineReIm j=%d beta=%d: %w", j, beta, err)
			}
			ci, err := c.Sw.StandardToCI(z)
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("refresh StandardToCI j=%d beta=%d: %w", j, beta, err)
			}
			ci.Scale = c.Canon
			out[j][beta] = ci
		}
		fmt.Printf("[refresh] recombine group %2d/8 done\n", j+1)
	}
	return out, nil
}

// AddRoundKey XORs the round key into the state, bit by bit. rk must already be at ARKLevel
func (c *Context) AddRoundKey(st, rk blockpack.Packed) (blockpack.Packed, error) {
	if kl, sl := rk[0][0].Level(), st[0][0].Level(); kl != sl {
		return blockpack.Packed{}, fmt.Errorf("AddRoundKey: round key at level %d, state at level %d", kl, sl)
	}
	var out blockpack.Packed
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			x, err := c.AE.Xor(st[g][b], rk[g][b])
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("AddRoundKey xor [%d][%d]: %w", g, b, err)
			}
			out[g][b] = x
		}
	}
	return out, nil
}

// Clean pulls every bit of the state back onto 0 and 1 after the round has spread them.
// Th second order cleaning function is used here(SmootherCleaning), see discussion in []
func (c *Context) Clean(st blockpack.Packed) (blockpack.Packed, error) {
	var out blockpack.Packed
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			cc, err := cleaning.SmootherCleaning(c.Sw.CiP, c.Sw.EvalCI, st[g][b])
			if err != nil {
				return blockpack.Packed{}, fmt.Errorf("Clean [%d][%d]: %w", g, b, err)
			}
			out[g][b] = cc
		}
	}
	return out, nil
}

// RoundStep is called after each operation of a middle round: its name, its duration, the state it
// produced. That is where a caller hooks its timings and its oracle checks; nil skips it.
type RoundStep func(name string, dur time.Duration, st blockpack.Packed)

// Round advances the packed state by one AES middle round: SubBytes -> Refresh -> MixColumns ->
// AddRoundKey -> Cleaning. rk is the round key at ARKLevel, version the SubByte variant (2..4).
func (c *Context) Round(st, rk blockpack.Packed, version int, after RoundStep) (blockpack.Packed, error) {
	if after == nil {
		after = func(string, time.Duration, blockpack.Packed) {}
	}
	var err error

	fmt.Println("---- SubBytes ----")
	t0 := time.Now()
	if st, err = c.SubBytes(st, version); err != nil {
		return blockpack.Packed{}, fmt.Errorf("Round SubBytes: %w", err)
	}
	after("SubBytes", time.Since(t0), st)

	fmt.Println("---- Refresh ----")
	t0 = time.Now()
	if st, err = c.Refresh(st); err != nil {
		return blockpack.Packed{}, fmt.Errorf("Round Refresh: %w", err)
	}
	after("Refresh", time.Since(t0), st)

	fmt.Println("---- MixColumns ----")
	t0 = time.Now()
	if st, err = c.AE.MixColumns(st); err != nil {
		return blockpack.Packed{}, fmt.Errorf("Round MixColumns: %w", err)
	}
	after("MixColumns", time.Since(t0), st)

	fmt.Println("---- AddRoundKey ----")
	t0 = time.Now()
	if st, err = c.AddRoundKey(st, rk); err != nil {
		return blockpack.Packed{}, fmt.Errorf("Round AddRoundKey: %w", err)
	}
	after("AddRoundKey", time.Since(t0), st)

	fmt.Println("---- Cleaning ----")
	t0 = time.Now()
	if st, err = c.Clean(st); err != nil {
		return blockpack.Packed{}, fmt.Errorf("Round Cleaning: %w", err)
	}
	after("Cleaning", time.Since(t0), st)

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

	fmt.Println("---- SubBytes ----")
	t0 := time.Now()
	if st, err = c.SubBytes(st, version); err != nil {
		return blockpack.Packed{}, fmt.Errorf("LastRoundV1 SubBytes: %w", err)
	}
	after("SubBytes", time.Since(t0), st)

	fmt.Println("---- Refresh ----")
	t0 = time.Now()
	if st, err = c.Refresh(st); err != nil {
		return blockpack.Packed{}, fmt.Errorf("LastRoundV1 Refresh: %w", err)
	}
	after("Refresh", time.Since(t0), st)

	fmt.Println("---- AddRoundKey ----")
	t0 = time.Now()
	if st, err = c.AddRoundKey(st, rk); err != nil {
		return blockpack.Packed{}, fmt.Errorf("LastRoundV1 AddRoundKey: %w", err)
	}
	after("AddRoundKey", time.Since(t0), st)

	fmt.Println("---- Cleaning ----")
	t0 = time.Now()
	if st, err = c.Clean(st); err != nil {
		return blockpack.Packed{}, fmt.Errorf("LastRoundV1 Cleaning: %w", err)
	}
	after("Cleaning", time.Since(t0), st)

	return st, nil
}
