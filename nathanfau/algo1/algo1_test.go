package algo1

import (
	"fmt"
	"math"
	"math/big"
	"math/cmplx"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/bitbatching"
	"github.com/tuneinsight/lattigo/v6/nathanfau/cleaning"
	"github.com/tuneinsight/lattigo/v6/nathanfau/convctx"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// newAlgo1Context builds the machinery to run Algo1. The parameters are built CI-first
// (ConjugateInvariant ring, then Standard) so the Q primes are NTT-friendly for the CI ring
// (= 1 mod 4N), which lets the SAME params serve both the slim bootstrap and the CI domain
// switch. One secret key is shared between the bootstrapping evaluator and the CtxSwitcher.
func newAlgo1Context(t *testing.T, k int) (ckks.Parameters, *ckks.Encoder, *rlwe.Encryptor, *rlwe.Decryptor, *bootstrapping.Evaluator, *convctx.CtxSwitcher) {
	t.Helper()

	logQ := []int{42}
	logQ = append(logQ, 60)                 // SlotsToCoeffs
	logQ = append(logQ, 38, 38)             // cleaning
	logQ = append(logQ, 38, 38, 38, 38, 38) // 5 levels for Bitextracting
	logQ = append(logQ, 38, 38, 38)         // 3 levels for squarring
	logQ = append(logQ, 38)                 // etractExp
	logQ = append(logQ, 38)                 // second Conv_{Real->Cplx}
	logQ = append(logQ, 38, 38, 38, 38, 38) // EvalCos		(38*7 in the article)
	logQ = append(logQ, 38)                 // first Conv_{Cplx->Real}
	logQ = append(logQ, 38, 38, 38)         // CoeffToSlots

	ciBase, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            11,
		LogQ:            logQ,
		LogP:            []int{61, 61, 61, 61, 61, 61, 50},
		LogDefaultScale: 38,
		Xs:              ring.Ternary{H: 256},
		RingType:        ring.ConjugateInvariant,
	})
	if err != nil {
		t.Fatalf("ci params: %v", err)
	}
	stdLit := ciBase.ParametersLiteral()
	stdLit.RingType = ring.Standard
	params, err := ckks.NewParametersFromLiteral(stdLit)
	if err != nil {
		t.Fatalf("std params: %v", err)
	}

	// SlotsToCoeffs (homomorphic decoding).
	S2CParams := dft.MatrixLiteral{
		Type:     dft.HomomorphicDecode,
		LogSlots: params.LogMaxSlots(),
		LevelP:   params.MaxLevelP(),
		Levels:   []int{1},
	}
	S2CParams.LevelQ = len(S2CParams.Levels)

	// CoeffsToSlots (homomorphic encoding), real and imaginary parts split.
	C2SParams := dft.MatrixLiteral{
		Type:     dft.HomomorphicEncode,
		Format:   dft.SplitRealAndImag,
		LogSlots: params.LogMaxSlots(),
		LevelQ:   params.MaxLevel(),
		LevelP:   params.MaxLevelP(),
		Levels:   []int{1, 1, 1},
	}

	// Mod1 parameters, mostly unused (EvalMod is replaced by the trig step) but still needed to
	// build the bootstrapping keys and expose ScalingFactor and MessageRatio.
	Mod1Params := mod1.ParametersLiteral{
		LevelQ:          params.MaxLevel() - C2SParams.Depth(true),
		LogScale:        38,
		Mod1Type:        mod1.CosDiscrete,
		Mod1Degree:      2 * ((1 << k) - 1),
		K:               1 << k,
		LogMessageRatio: k,
	}
	mod1P, err := mod1.NewParametersFromLiteral(params, Mod1Params)
	if err != nil {
		t.Fatalf("mod1 params: %v", err)
	}
	F := (mod1P.ScalingFactor().Float64() / mod1P.MessageRatio()) / params.DefaultScale().Float64()
	S2CParams.Scaling = big.NewFloat(F)

	btpParams := bootstrapping.Parameters{
		ResidualParameters:      params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: S2CParams,
		Mod1ParametersLiteral:   Mod1Params,
		CoeffsToSlotsParameters: C2SParams,
		EphemeralSecretWeight:   32,
		CircuitOrder:            bootstrapping.DecodeThenModUp,
	}

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

	ecd := ckks.NewEncoder(params)
	enc := rlwe.NewEncryptor(params, sk)
	dec := rlwe.NewDecryptor(params, sk)

	debug.EncStd, debug.DecStd, debug.ParamsStd = ecd, dec, params
	debug.EncCI, debug.DecCI, debug.ParamsCI = ckks.NewEncoder(sw.CiP), rlwe.NewDecryptor(sw.CiP, sw.SkCI), sw.CiP

	return params, ecd, enc, dec, eval, sw
}

// TestAlgo1 runs the bootstrap end to end (no transciphering, no AES). It encrypts k
// bit-planes of two independent streams packed as the real and imaginary parts, packs them into
// one integer per stream with BitPack, bootstraps with Algo1 (which returns the t-th roots of
// unity of each stream), extracts the bits back with BitExtract, cleans them, and checks that
// every recovered bit matches the original plane.
func TestAlgo1(t *testing.T) {
	const k = 4
	params, ecd, enc, dec, eval, sw := newAlgo1Context(t, k)
	nSlots := params.MaxSlots()
	tMod := 1 << k
	level := eval.SlotsToCoeffsParameters.LevelQ

	rng := rand.New(rand.NewSource(time.Now().UnixNano())) // *random* seed for testing

	// k bit-planes, two independent streams packed as the real and imaginary parts.
	bitsRe := make([][]float64, k)
	bitsIm := make([][]float64, k)
	have := make([]*rlwe.Ciphertext, k)
	for j := 0; j < k; j++ {
		bitsRe[j] = make([]float64, nSlots)
		bitsIm[j] = make([]float64, nSlots)
		vals := make([]complex128, nSlots)
		for s := 0; s < nSlots; s++ {
			bitsRe[j][s] = float64(rng.Intn(2))
			bitsIm[j][s] = float64(rng.Intn(2))
			vals[s] = complex(bitsRe[j][s], bitsIm[j][s])
		}
		pt := ckks.NewPlaintext(params, level)
		if err := ecd.Encode(vals, pt); err != nil {
			t.Fatalf("encode bit-plane %d: %v", j, err)
		}
		ct, err := enc.EncryptNew(pt)
		if err != nil {
			t.Fatalf("encrypt bit-plane %d: %v", j, err)
		}
		have[j] = ct
	}

	fmt.Println("================ input bit-planes ================")
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
	fmt.Println("================ BitPack ================")
	debug.DbgSlotStd("ctPack =", ctPack)
	debug.DbgChain("chain ctPack :", eval.Evaluator, ctPack)

	// Expected roots of unity per stream (precision traces only).
	wantRe := make([]complex128, nSlots)
	wantIm := make([]complex128, nSlots)
	for s := 0; s < nSlots; s++ {
		var mRe, mIm float64
		for j := 0; j < k; j++ {
			mRe += bitsRe[j][s] * math.Exp2(float64(j))
			mIm += bitsIm[j][s] * math.Exp2(float64(j))
		}
		wantRe[s] = cmplx.Exp(complex(0, 2*math.Pi*mRe/float64(tMod)))
		wantIm[s] = cmplx.Exp(complex(0, 2*math.Pi*mIm/float64(tMod)))
	}

	// Algo1 = the bootstrap.
	fmt.Println("================ Algo1 (bootstrap) ================")
	ctReal, ctImag, err := Algo1(eval, sw, ctPack, k)
	if err != nil {
		t.Fatalf("Algo1: %v", err)
	}
	debug.DbgSlotStd("ctReal =", ctReal)
	debug.DbgChain("chain ctReal :", eval.Evaluator, ctReal)
	debug.PrintPrecStd("prec ctReal (roots) :", wantRe, ctReal)
	debug.DbgSlotStd("ctImag =", ctImag)
	debug.DbgChain("chain ctImag :", eval.Evaluator, ctImag)
	debug.PrintPrecStd("prec ctImag (roots) :", wantIm, ctImag)

	// BitExtract each stream back to its k bits.
	bitsReExtrated, err := bitbatching.BitExtract(params, eval.Evaluator, ctReal, k)
	if err != nil {
		t.Fatalf("BitExtract real: %v", err)
	}
	bitsImExtrated, err := bitbatching.BitExtract(params, eval.Evaluator, ctImag, k)
	if err != nil {
		t.Fatalf("BitExtract imag: %v", err)
	}
	fmt.Println("================ BitExtract (raw) ================")
	for j := 0; j < k; j++ {
		debug.DbgSlotStd(fmt.Sprintf("bitsReExtrated[%d] raw =", j), bitsReExtrated[j])
		debug.DbgChain("chain bitsReExtrated raw :", eval.Evaluator, bitsReExtrated[j])
		debug.PrintPrecStd("prec bitsReExtrated raw :", bitsRe[j], bitsReExtrated[j])
		debug.DbgSlotStd(fmt.Sprintf("bitsImExtrated[%d] raw =", j), bitsImExtrated[j])
		debug.DbgChain("chain bitsImExtrated raw :", eval.Evaluator, bitsImExtrated[j])
		debug.PrintPrecStd("prec bitsImExtrated raw :", bitsIm[j], bitsImExtrated[j])
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
	fmt.Println("================ Cleaning ================")
	for j := 0; j < k; j++ {
		debug.DbgSlotStd(fmt.Sprintf("bitsReExtrated[%d] cleaned =", j), bitsReExtrated[j])
		debug.DbgChain("chain bitsReExtrated cleaned :", eval.Evaluator, bitsReExtrated[j])
		debug.PrintPrecStd("prec bitsReExtrated cleaned :", bitsRe[j], bitsReExtrated[j])
		debug.DbgSlotStd(fmt.Sprintf("bitsImExtrated[%d] cleaned =", j), bitsImExtrated[j])
		debug.DbgChain("chain bitsImExtrated cleaned :", eval.Evaluator, bitsImExtrated[j])
		debug.PrintPrecStd("prec bitsImExtrated cleaned :", bitsIm[j], bitsImExtrated[j])
	}

	// Check every recovered bit against the original plane.
	checkBits(t, ecd, dec, bitsReExtrated, bitsRe, "real")
	checkBits(t, ecd, dec, bitsImExtrated, bitsIm, "imag")
}

// checkBits decodes each extracted bit ciphertext and checks it against the original plane
func checkBits(t *testing.T, ecd *ckks.Encoder, dec *rlwe.Decryptor, got []*rlwe.Ciphertext, want [][]float64, name string) {
	t.Helper()
	for j := range got {
		dcd := make([]complex128, got[j].Slots())
		if err := ecd.Decode(dec.DecryptNew(got[j]), dcd); err != nil {
			t.Fatalf("%s decode bit %d: %v", name, j, err)
		}
		wrong := 0
		for s := range want[j] {
			if (real(dcd[s]) >= 0.5) != (want[j][s] >= 0.5) {
				wrong++
			}
		}
		t.Logf("%s bit %d: %d/%d wrong", name, j, wrong, len(want[j]))
		if wrong != 0 {
			t.Errorf("%s bit %d: %d bits wrong", name, j, wrong)
		}
	}
}
