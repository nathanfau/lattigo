package algo1

//	go test ./nathanfau/algo1/ -run '^TestAlgo1Zone$' -v -zone 33 -timeout 0
//	go test ./nathanfau/algo1/ -run '^TestAlgo1Zone$' -v -zone 33 -pre 3 -post 4 -timeout 0

import (
	"flag"
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
	"github.com/tuneinsight/lattigo/v6/nathanfau/bitbatching"
	"github.com/tuneinsight/lattigo/v6/nathanfau/cleaning"
	"github.com/tuneinsight/lattigo/v6/nathanfau/convctx"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/nathanfau/params2"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var (
	zoneFlag = flag.Int("zone", 33, "size of the zone's primes and scale (params2.Shape.LogZone); 38 is the chain without a zone")
	preFlag  = flag.Int("pre", 0, "XOR layers before the refresh, one level each (3 = SubBytes' depth); 0 = none, as in TestAlgo1")
	postFlag = flag.Int("post", 0, "XOR layers between BitExtract and the cleaning, one level each (4 = MixColumns + AddRoundKey); 0 = none, as in TestAlgo1")
)

const (
	// zoneScaleTol is how close, in bits, a ciphertext must sit to where an edge of the zone has to
	// leave it. Off by more means its scale was rounded somewhere, i.e. its value is wrong.
	zoneScaleTol = 30
	// zoneBitPrec is the bar the recovered bits must clear: |got - bit| <= 2^-zoneBitPrec. It asks
	// that the bits come back, not that they come back precise (TestAlgo1 asks 2^-20); that is for
	// later.
	zoneBitPrec = 10.0
)

// zoneAt checks that got sits at want, the scale an edge of the zone must leave the ciphertext at.
func zoneAt(t *testing.T, where string, got, want rlwe.Scale) {
	t.Helper()
	d := got.Log2Delta(want)
	fmt.Printf("  %-40s scale 2^%.6f, want 2^%.6f (off by 2^-%.1f)\n", where, got.Log2(), want.Log2(), d)
	if d < zoneScaleTol {
		t.Fatalf("%s: scale 2^%.6f is 2^-%.1f off 2^%.6f, want below 2^-%d", where, got.Log2(), d, want.Log2(), zoneScaleTol)
	}
}

// newZoneContext builds the machinery on the transciphering chain itself, with its zone at logZone
// bits, the way transciphering.NewContext does.
func newZoneContext(t *testing.T, logN, k, logZone int) (ckks.Parameters, rlwe.Scale, *bootstrapping.Evaluator, *convctx.CtxSwitcher, *rlwe.SecretKey) {
	t.Helper()

	sh := params2.Shape{LogZone: logZone}
	params, btpParams, err := params2.TranscipheringParamsWith(logN, k, 2, k, sh)
	if err != nil {
		t.Fatalf("TranscipheringParamsWith: %v", err)
	}
	debug.DbgParams(fmt.Sprintf("algo1 zone %d", logZone), params)

	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		t.Fatalf("GenEvaluationKeys: %v", err)
	}
	eval, err := bootstrapping.NewEvaluator(btpParams, evk)
	if err != nil {
		t.Fatalf("bootstrapping NewEvaluator: %v", err)
	}
	sw, err := convctx.NewCtxSwitcher(params, sk)
	if err != nil {
		t.Fatalf("NewCtxSwitcher: %v", err)
	}
	return params, sh.ZoneScale(), eval, sw, sk
}

// zoneKit is what the helpers of TestAlgo1Zone share: encryption in Std at the zone's scale, the
// pipeline's XOR, and the draw of fresh bits.
type zoneKit struct {
	params ckks.Parameters
	ecd    *ckks.Encoder
	enc    *rlwe.Encryptor
	xe     *aes.Evaluator
	rng    *rand.Rand
	zone   rlwe.Scale
}

// bits draws n uniform bits.
func (z *zoneKit) bits(n int) []float64 {
	b := make([]float64, n)
	for s := range b {
		b[s] = float64(z.rng.Intn(2))
	}
	return b
}

// encrypt encrypts vals ([]float64 or []complex128) at the given level, at the zone's scale.
func (z *zoneKit) encrypt(t *testing.T, vals any, level int) *rlwe.Ciphertext {
	t.Helper()
	pt := ckks.NewPlaintext(z.params, level)
	pt.Scale = z.zone
	if err := z.ecd.Encode(vals, pt); err != nil {
		t.Fatalf("encode at level %d: %v", level, err)
	}
	out, err := z.enc.EncryptNew(pt)
	if err != nil {
		t.Fatalf("encrypt at level %d: %v", level, err)
	}
	return out
}

// xorLayers stands in for the AES circuit around the refresh: n layers of aes.XorNoSq, one level
// each, on the k real bit-planes of one stream, want following in the clear.
func (z *zoneKit) xorLayers(t *testing.T, name string, n int, cts []*rlwe.Ciphertext, want [][]float64) {
	t.Helper()
	k := len(cts)
	for i := 0; i < n; i++ {
		next := make([]*rlwe.Ciphertext, k)
		nextWant := make([][]float64, k)
		for j := 0; j < k; j++ {
			partner, pw := cts[(j+1)%k], want[(j+1)%k]
			if (i+j)%2 == 0 {
				pw = z.bits(len(want[j]))
				partner = z.encrypt(t, pw, cts[j].Level())
			}
			x, err := z.xe.XorNoSq(cts[j], partner)
			if err != nil {
				t.Fatalf("%s layer %d plane %d: %v", name, i+1, j, err)
			}
			next[j] = x
			nextWant[j] = make([]float64, len(pw))
			for s := range pw {
				nextWant[j][s] = float64(int(want[j][s]) ^ int(pw[s]))
			}
		}
		copy(cts, next)
		copy(want, nextWant)

		line := fmt.Sprintf("  [%s %d] lv %2d  scale/zone - 1 :", name, i+1, cts[0].Level())
		for _, ct := range cts {
			line += fmt.Sprintf(" %+.3e", ct.Scale.Div(z.zone).Float64()-1)
		}
		fmt.Println(line)
	}
	spread := math.Inf(1)
	for _, ct := range cts[1:] {
		spread = math.Min(spread, ct.Scale.Log2Delta(cts[0].Scale))
	}
	fmt.Printf("  [%s] the %d planes' scales differ by up to 2^-%.1f\n", name, k, spread)
}

// TestAlgo1Zone is TestAlgo1, step for step and print for print, on the transciphering chain with
// its zone at -zone bits (params2.Shape.LogZone). The bit-planes are encrypted at the zone's scale;
// they leave the zone for the bootstrap through a scale hop onto EntryScale, on the
// Conv_{Real->Cplx} prime (the one CIToStandardTo spends in the pipeline); ExtractTo brings
// Algo1's output back onto the zone's scale, where BitExtract and the cleaning run on the zone's
// primes. The scale is checked at each edge (the "scale ... want ..." lines, the only ones TestAlgo1
// does not print).
//
// -pre and -post add XOR layers before the refresh and after BitExtract, standing in for the AES
// circuit: they move the scale off 2^LogZone the way the real circuit does, which the edges have to
// absorb. XORs are real, so with -pre the two streams are encrypted and XORed apart, then
// recombined as Re + i.Im; the post layers run on BitExtract's outputs, which are already apart.
func TestAlgo1Zone(t *testing.T) {
	const logN, k = 11, 4
	params, zone, eval, sw, sk := newZoneContext(t, logN, k, *zoneFlag)
	ecd := ckks.NewEncoder(params)
	dec := rlwe.NewDecryptor(params, sk)
	debug.EncStd, debug.DecStd, debug.ParamsStd = ecd, dec, params
	debug.EncCI, debug.DecCI, debug.ParamsCI = ckks.NewEncoder(sw.CiP), rlwe.NewDecryptor(sw.CiP, sw.SkCI), sw.CiP
	z := &zoneKit{
		params: params,
		ecd:    ecd,
		enc:    rlwe.NewEncryptor(params, sk),
		xe:     aes.NewEvaluator(eval.Evaluator),
		rng:    rand.New(rand.NewSource(time.Now().UnixNano())), // *random* seed for testing
		zone:   zone,
	}

	nSlots := params.MaxSlots()
	tMod := 1 << k
	// One prime above SlotsToCoeffs, the one the hop out of the zone rescales by, plus one per XOR.
	level := eval.SlotsToCoeffsParameters.LevelQ + 1 + *preFlag
	if level > params.MaxLevel() {
		t.Fatalf("-pre %d needs the input at level %d, the chain stops at %d", *preFlag, level, params.MaxLevel())
	}

	// k bit-planes, two independent streams packed as the real and imaginary parts.
	bitsRe := make([][]float64, k)
	bitsIm := make([][]float64, k)
	for j := 0; j < k; j++ {
		bitsRe[j] = z.bits(nSlots)
		bitsIm[j] = z.bits(nSlots)
	}
	have := make([]*rlwe.Ciphertext, k)
	if *preFlag == 0 {
		for j := 0; j < k; j++ {
			vals := make([]complex128, nSlots)
			for s := 0; s < nSlots; s++ {
				vals[s] = complex(bitsRe[j][s], bitsIm[j][s])
			}
			have[j] = z.encrypt(t, vals, level)
		}
	} else {
		haveRe := make([]*rlwe.Ciphertext, k)
		haveIm := make([]*rlwe.Ciphertext, k)
		for j := 0; j < k; j++ {
			haveRe[j] = z.encrypt(t, bitsRe[j], level)
			haveIm[j] = z.encrypt(t, bitsIm[j], level)
		}
		fmt.Println("---- XOR pre ----")
		z.xorLayers(t, "pre Re", *preFlag, haveRe, bitsRe)
		z.xorLayers(t, "pre Im", *preFlag, haveIm, bitsIm)
		for j := 0; j < k; j++ {
			var err error
			if have[j], err = utils.CombineReIm(eval.Evaluator, haveRe[j], haveIm[j]); err != nil {
				t.Fatalf("CombineReIm plane %d: %v", j, err)
			}
		}
	}

	fmt.Println("----input bit-planes----")
	for j := 0; j < k; j++ {
		debug.DbgSlotStd(fmt.Sprintf("have[%d] =", j), have[j])
		debug.DbgChain("chain have :", eval.Evaluator, have[j])
		debug.PrintPrecStd("prec have :", bitsRe[j], have[j])
	}

	// BitPack the k planes into one packed integer per stream, complex(m_re, m_im).
	groups, err := bitbatching.BitPackGroups(eval.Evaluator, have, k)
	if err != nil {
		t.Fatalf("BitPackGroups: %v", err)
	}
	if len(groups) != 1 {
		t.Fatalf("expected 1 group, got %d", len(groups))
	}
	ctPack := groups[0]
	fmt.Println("---- BitPack ----")
	debug.DbgSlotStd("ctPack =", ctPack)
	debug.DbgChain("chain ctPack :", eval.Evaluator, ctPack)

	// Expected packed integers and roots of unity per stream (precision traces only).
	wantPack := make([]complex128, nSlots)
	wantRe := make([]complex128, nSlots)
	wantIm := make([]complex128, nSlots)
	for s := 0; s < nSlots; s++ {
		var mRe, mIm float64
		for j := 0; j < k; j++ {
			mRe += bitsRe[j][s] * math.Exp2(float64(j))
			mIm += bitsIm[j][s] * math.Exp2(float64(j))
		}
		wantPack[s] = complex(mRe, mIm)
		wantRe[s] = cmplx.Exp(complex(0, 2*math.Pi*mRe/float64(tMod)))
		wantIm[s] = cmplx.Exp(complex(0, 2*math.Pi*mIm/float64(tMod)))
	}
	debug.PrintPrecStd("prec ctPack :", wantPack, ctPack)

	// Out of the zone: the scale hop onto what Algo1 wants, on the Conv_{Real->Cplx} prime.
	entry := EntryScale(eval)
	if err = utils.MulRescaleTo(eval.Evaluator, ctPack, 1, entry); err != nil {
		t.Fatalf("hop out of the zone: %v", err)
	}
	zoneAt(t, "hop out of the zone", ctPack.Scale, entry)

	// Algo1 = the bootstrap, landing back on the zone's scale.
	fmt.Println("---- Algo1 (bootstrap) ----")
	ctReal, ctImag, err := ExtractTo(eval, sw, ctPack, k, zone)
	if err != nil {
		t.Fatalf("ExtractTo: %v", err)
	}
	if err = Resume(eval, ctReal, ctImag); err != nil {
		t.Fatalf("Resume: %v", err)
	}
	debug.DbgSlotStd("ctReal =", ctReal)
	debug.DbgChain("chain ctReal :", eval.Evaluator, ctReal)
	debug.PrintPrecStd("prec ctReal (roots) :", wantRe, ctReal)
	debug.DbgSlotStd("ctImag =", ctImag)
	debug.DbgChain("chain ctImag :", eval.Evaluator, ctImag)
	debug.PrintPrecStd("prec ctImag (roots) :", wantIm, ctImag)
	zoneAt(t, "Algo1 real (back in the zone)", ctReal.Scale, zone)
	zoneAt(t, "Algo1 imag (back in the zone)", ctImag.Scale, zone)

	// BitExtract each stream back to its k bits.
	bitsReExtrated, err := bitbatching.BitExtract(params, eval.Evaluator, ctReal, k)
	if err != nil {
		t.Fatalf("BitExtract real: %v", err)
	}
	bitsImExtrated, err := bitbatching.BitExtract(params, eval.Evaluator, ctImag, k)
	if err != nil {
		t.Fatalf("BitExtract imag: %v", err)
	}
	fmt.Println("---- BitExtract  ----")
	for j := 0; j < k; j++ {
		debug.DbgSlotStd(fmt.Sprintf("bitsReExtrated[%d] raw =", j), bitsReExtrated[j])
		debug.DbgChain("chain bitsReExtrated raw :", eval.Evaluator, bitsReExtrated[j])
		debug.PrintPrecStd("prec bitsReExtrated raw :", bitsRe[j], bitsReExtrated[j])
		debug.DbgSlotStd(fmt.Sprintf("bitsImExtrated[%d] raw =", j), bitsImExtrated[j])
		debug.DbgChain("chain bitsImExtrated raw :", eval.Evaluator, bitsImExtrated[j])
		debug.PrintPrecStd("prec bitsImExtrated raw :", bitsIm[j], bitsImExtrated[j])
	}

	if *postFlag > 0 {
		if lv := bitsReExtrated[0].Level() - *postFlag; lv < 2 {
			t.Fatalf("-post %d leaves level %d, the cleaning needs 2", *postFlag, lv)
		}
		fmt.Println("---- XOR post ----")
		z.xorLayers(t, "post Re", *postFlag, bitsReExtrated, bitsRe)
		z.xorLayers(t, "post Im", *postFlag, bitsImExtrated, bitsIm)
		for j := 0; j < k; j++ {
			debug.PrintPrecStd("prec bitsReExtrated post :", bitsRe[j], bitsReExtrated[j])
			debug.PrintPrecStd("prec bitsImExtrated post :", bitsIm[j], bitsImExtrated[j])
		}
	}

	// Cleaning snaps the extracted bits to 0/1.
	for j := 0; j < k; j++ {
		if bitsReExtrated[j], err = cleaning.Cleaning(params, eval.Evaluator, bitsReExtrated[j]); err != nil {
			t.Fatalf("Cleaning real bit %d: %v", j, err)
		}
		if bitsImExtrated[j], err = cleaning.Cleaning(params, eval.Evaluator, bitsImExtrated[j]); err != nil {
			t.Fatalf("Cleaning imag bit %d: %v", j, err)
		}
	}
	fmt.Println("---- Cleaning ----")
	for j := 0; j < k; j++ {
		debug.DbgSlotStd(fmt.Sprintf("bitsReExtrated[%d] cleaned =", j), bitsReExtrated[j])
		debug.DbgChain("chain bitsReExtrated cleaned :", eval.Evaluator, bitsReExtrated[j])
		debug.PrintPrecStd("prec bitsReExtrated cleaned :", bitsRe[j], bitsReExtrated[j])
		debug.DbgSlotStd(fmt.Sprintf("bitsImExtrated[%d] cleaned =", j), bitsImExtrated[j])
		debug.DbgChain("chain bitsImExtrated cleaned :", eval.Evaluator, bitsImExtrated[j])
		debug.PrintPrecStd("prec bitsImExtrated cleaned :", bitsIm[j], bitsImExtrated[j])
	}

	// Check every recovered bit against the original plane.
	check := func(name string, got []*rlwe.Ciphertext, want [][]float64) {
		res, err := utils.BitDistance(ecd, dec, got, want, math.Exp2(-zoneBitPrec))
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		for j, r := range res {
			msg := fmt.Sprintf("%s bit %d: worst |err| = 2^%.1f, %d/%d slots over 2^%.0f",
				name, j, math.Log2(r.Worst), r.Wrong, r.Slots, math.Log2(r.Tol))
			t.Log(msg)
			if r.Wrong != 0 {
				t.Error(msg)
			}
		}
	}
	check("real", bitsReExtrated, bitsRe)
	check("imag", bitsImExtrated, bitsIm)
}
