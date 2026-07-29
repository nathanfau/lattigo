// Package params holds the CKKS and bootstrapping parameter sets for the transciphering tests
// (nathanfau/transciphering): the k=4 Algo1 pipeline (TranscipheringParams). The other packages'
// tests build their own parameters.
package params2

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// ciThenStd builds a conjugate-invariant parameter set from the given moduli chain, then its
// standard-ring twin (same Q, P). The CI-first order makes the Q primes NTT-friendly for the CI
// ring (= 1 mod 4N), so the same primes serve both the bootstrap and the CI domain switch.
func ciThenStd(logN int, logQ, logP []int, logScale int) (ckks.Parameters, error) {
	ciBase, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            logN,
		LogQ:            logQ,
		LogP:            logP,
		LogDefaultScale: logScale,
		Xs:              ring.Ternary{H: 256},
		RingType:        ring.ConjugateInvariant,
	})
	if err != nil {
		return ckks.Parameters{}, fmt.Errorf("ci params: %w", err)
	}
	stdLit := ciBase.ParametersLiteral()
	stdLit.RingType = ring.Standard
	params, err := ckks.NewParametersFromLiteral(stdLit)
	if err != nil {
		return ckks.Parameters{}, fmt.Errorf("std params: %w", err)
	}
	return params, nil
}

// TranscipheringParams is the Algo1 pipeline parameter set (the leveled AES round + Algo1 refresh
// of nathanfau/transciphering): MessageRatio = 2^k so q0 = base + k = 42 at k = 4, K = t = 2^k.
func TranscipheringParams(logN, k int) (ckks.Parameters, bootstrapping.Parameters, error) {
	logQ := []int{42}
	logQ = append(logQ, 60)         // SlotsToCoeffs
	logQ = append(logQ, 38)         // Conv_{Real->Cplx}
	logQ = append(logQ, 38, 38, 38) // SubBytes
	logQ = append(logQ, 38, 38)     // cleaning
	logQ = append(logQ, 38)         // + 1 additional level to run SmootherCleaning or VerySmoothrCleaning
	logQ = append(logQ, 38)         // AddRoundKey
	logQ = append(logQ, 38, 38, 38) // MixColumns
	// refresh
	logQ = append(logQ, 38)                 // Conv_{Cplx->Real}
	logQ = append(logQ, 38, 38, 38, 38, 38) // 5 levels for BitExtract
	logQ = append(logQ, 38, 38, 38)         // 3 levels for squaring
	logQ = append(logQ, 38)                 // extractExp
	logQ = append(logQ, 38)                 // Conv_{Real->Cplx}
	logQ = append(logQ, 38, 38, 38, 38, 38) // EvalCos
	logQ = append(logQ, 38)                 // Conv_{Cplx->Real}
	logQ = append(logQ, 38, 38, 38)         // CoeffsToSlots

	params, err := ciThenStd(logN, logQ, []int{61, 61, 61, 61, 61, 61, 50}, 38)
	if err != nil {
		return ckks.Parameters{}, bootstrapping.Parameters{}, err
	}

	S2CParams := dft.MatrixLiteral{
		Type:     dft.HomomorphicDecode,
		LogSlots: params.LogMaxSlots(),
		LevelP:   params.MaxLevelP(),
		Levels:   []int{1},
	}
	S2CParams.LevelQ = len(S2CParams.Levels)

	C2SParams := dft.MatrixLiteral{
		Type:     dft.HomomorphicEncode,
		Format:   dft.SplitRealAndImag,
		LogSlots: params.LogMaxSlots(),
		LevelQ:   params.MaxLevel(),
		LevelP:   params.MaxLevelP(),
		Levels:   []int{1, 1, 1},
	}

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
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("mod1 params: %w", err)
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
	return params, btpParams, nil
}
