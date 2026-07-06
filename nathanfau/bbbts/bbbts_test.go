package bbbts

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
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// newBootstrapContext builds a slim-order bootstrapping evaluator for k-bit packing.
func newBootstrapContext(t *testing.T, k int) (ckks.Parameters, *ckks.Encoder, *rlwe.Encryptor, *rlwe.Decryptor, *bootstrapping.Evaluator) {
	t.Helper()

	logQ := []int{50}
	for i := 0; i < 19; i++ {
		logQ = append(logQ, 45)
	}

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            11,
		LogQ:            logQ,
		LogP:            []int{61, 61, 61, 61, 61, 61, 50},
		LogDefaultScale: 45,
		Xs:              ring.Ternary{H: 256},
	})
	if err != nil {
		t.Fatalf("params: %v", err)
	}

	// SlotsToCoeffs (homomorphic decoding).
	S2CParams := dft.MatrixLiteral{
		Type:     dft.HomomorphicDecode,
		LogSlots: params.LogMaxSlots(),
		LevelP:   params.MaxLevelP(),
		Levels:   []int{1},
	}
	S2CParams.LevelQ = len(S2CParams.Levels)

	// CoeffsToSlots (homomorphic encoding); splits real and imaginary parts.
	C2SParams := dft.MatrixLiteral{
		Type:     dft.HomomorphicEncode,
		Format:   dft.SplitRealAndImag,
		LogSlots: params.LogMaxSlots(),
		LevelQ:   params.MaxLevel(),
		LevelP:   params.MaxLevelP(),
		Levels:   []int{1, 1, 1},
	}

	// Mod1 parameters: mostly unused (EvalMod is replaced by EvalExp) but still needed to
	// build the bootstrapping keys and to expose ScalingFactor / MessageRatio.
	Mod1Params := mod1.ParametersLiteral{
		LevelQ:          params.MaxLevel() - C2SParams.Depth(true),
		LogScale:        45,
		Mod1Type:        mod1.CosDiscrete,
		Mod1Degree:      2 * ((1 << k) - 1),
		K:               1 << k,
		LogMessageRatio: k,
	}

	mod1P, err := mod1.NewParametersFromLiteral(params, Mod1Params)
	if err != nil {
		t.Fatalf("mod1 params: %v", err)
	}

	// Slim order: cancel the SlotsToCoeffs rescaling that normally follows EvalMod.
	// initialize() multiplies this by scale/offset = 1/F, so the net StC scaling is 1.
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

	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		t.Fatalf("GenEvaluationKeys: %v", err)
	}

	ecd := ckks.NewEncoder(params)
	enc := rlwe.NewEncryptor(params, pk)
	dec := rlwe.NewDecryptor(params, sk)

	debug.EncStd, debug.DecStd, debug.ParamsStd = ecd, dec, params

	eval, err := bootstrapping.NewEvaluator(btpParams, evk)
	if err != nil {
		t.Fatalf("NewEvaluator: %v", err)
	}

	return params, ecd, enc, dec, eval
}

// TestBBBTS runs the full BBBTS algorithm from [BKSS24]
func TestBBBTS(t *testing.T) {
	const k = 4
	const nGroups = 1
	nCt := nGroups * k

	params, ecd, enc, dec, eval := newBootstrapContext(t, k)
	nSlots := params.MaxSlots()
	level := eval.SlotsToCoeffsParameters.LevelQ

	rng := rand.New(rand.NewSource(time.Now().UnixNano())) // *random* seed for testing
	bits := make([][]float64, nCt)
	have := make([]*rlwe.Ciphertext, nCt)
	for iter := 0; iter < nCt; iter++ {
		b := make([]float64, nSlots)
		for s := range b {
			b[s] = float64(rng.Intn(2))
		}
		bits[iter] = b

		pt := ckks.NewPlaintext(params, level)
		if err := ecd.Encode(b, pt); err != nil {
			t.Fatalf("encode ct %d: %v", iter, err)
		}
		ct, err := enc.EncryptNew(pt)
		if err != nil {
			t.Fatalf("encrypt ct %d: %v", iter, err)
		}
		have[iter] = ct

		fmt.Printf("bits_%d: %f ...\n", iter, b[:2])
	}

	for iter := 0; iter < nCt; iter++ {
		debug.DbgSlotStd("ct_bits", have[iter])
		debug.DbgChain("chain ct_bits :", eval.Evaluator, have[iter])
	}

	for iter := 0; iter < nCt; iter++ {
		debug.PrintPrecStd("prec ct[i] :", bits[iter], have[iter])
	}

	// Bit packing by groups of k.
	groups, err := bitbatching.BitPackGroups(eval.Evaluator, have, k)
	if err != nil {
		t.Fatalf("BitPackGroups: %v", err)
	}

	tMod := 1 << k
	for g, ctPack := range groups {
		// Expected packed integer m = sum_i bit_i * 2^i, and its root of unity exp(2*pi*i*m/t).
		packed := make([]float64, nSlots)
		expWant := make([]complex128, nSlots)
		for s := 0; s < nSlots; s++ {
			for i := 0; i < k; i++ {
				packed[s] += bits[g*k+i][s] * math.Exp2(float64(i))
			}
			expWant[s] = cmplx.Exp(complex(0, 2*math.Pi*packed[s]/float64(tMod)))
		}

		debug.DbgSlotStd("ctPack =", ctPack)
		debug.DbgChain("chain ctPack :", eval.Evaluator, ctPack)
		debug.PrintPrecStd("prec ct_pack :", packed, ctPack)

		// IntRootBoot : slim bootstrapping with EvalExp instead of EvalMod.
		ctExp, err := IntRootBoot(eval, params, ctPack, k)
		if err != nil {
			t.Fatalf("group %d IntRootBoot: %v", g, err)
		}
		debug.DbgSlotStd("ctExp =", ctExp)
		debug.DbgChain("chain ctExp :", eval.Evaluator, ctExp)
		debug.PrintPrecStd("prec ctExp :", expWant, ctExp)

		// BitExtract
		ctBits, err := bitbatching.BitExtract(params, eval.Evaluator, ctExp, k)
		if err != nil {
			t.Fatalf("group %d BitExtract: %v", g, err)
		}

		for iter := 0; iter < k; iter++ {
			debug.DbgSlotStd("ct_bits (raw) =", ctBits[iter])
			debug.DbgChain("chain ct_bits (raw):", eval.Evaluator, ctBits[iter])
		}

		for iter := 0; iter < k; iter++ {
			debug.PrintPrecStd("prec bit (raw) :", bits[g*k+iter], ctBits[iter])
		}

		// Cleaning
		for iter := 0; iter < k; iter++ {
			ctBits[iter], err = cleaning.Cleaning(params, eval.Evaluator, ctBits[iter])
			if err != nil {
				t.Fatalf("group %d Cleaning bit %d: %v", g, iter, err)
			}
		}

		for iter := 0; iter < k; iter++ {
			debug.DbgSlotStd("ct_bits (cleaned) =", ctBits[iter])
			debug.DbgChain("chain ct_bits (cleaned) :", eval.Evaluator, ctBits[iter])
		}

		for iter := 0; iter < k; iter++ {
			debug.PrintPrecStd("prec bit (cleaned) :", bits[g*k+iter], ctBits[iter])
		}

		// Check each recovered bit against the original bit-plane bits[g*k+iter].
		for iter := 0; iter < k; iter++ {
			want := bits[g*k+iter]

			got := make([]complex128, ctBits[iter].Slots())
			if err := ecd.Decode(dec.DecryptNew(ctBits[iter]), got); err != nil {
				t.Fatalf("group %d decode bit %d: %v", g, iter, err)
			}
			maxErr := 0.0
			for s := 0; s < nSlots; s++ {
				if e := math.Abs(real(got[s]) - want[s]); e > maxErr {
					maxErr = e
				}
			}
			t.Logf("group %d bit %d: max err=%.3e", g, iter, maxErr)
			if maxErr > 1e-2 {
				t.Fatalf("group %d bit %d error %.3e too large", g, iter, maxErr)
			}
		}
	}
}
