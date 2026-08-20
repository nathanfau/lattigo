package params2

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// The parameter set of the Walsh-Hadamard pipeline (nathanfau/whtransciphering): a leveled AES
// round refreshed by SlotsToCoeffs, ModUp, CoeffsToSlots and a homomorphic x mod 1, in that order.
const (
	// WHLogScale is the working scale AND the size of the round primes. The XOR tree and
	// AddRoundKey are ct x ct, whose scale follows s -> 2s - p, and p = s is its only fixed point.
	WHLogScale = 28

	WHLogMessageRatio = 6

	// WHLogQ0 is three things at once, which is why it is written once: the size of q0, the scale
	// x mod 1 is held at, and the size of the EvalMod primes.
	WHLogQ0 = WHLogScale + WHLogMessageRatio + 1
)

type WHChain struct {
	Clean bool // two more primes, for a (3x-x^3)/2 on the round output

	C2SDepth    int // CoeffsToSlots levels
	S2CDepth    int // SlotsToCoeffs levels
	Mod1Degree  int // degree of the x mod 1 polynomial
	DoubleAngle int // double-angle steps of EvalMod
	K           int // half-width of the x mod 1 approximation interval
}

// DefaultWHChain is the set the pipeline runs on.
func DefaultWHChain() WHChain {
	return WHChain{
		C2SDepth:    4,
		S2CDepth:    2,
		Mod1Degree:  31,
		DoubleAngle: 3,
		K:           16,
	}
}

// WHTranscipheringParams builds the single moduli chain the whole pipeline runs on. The circuit
// order is DecodeThenModUp, so the chain is consumed once, top to bottom, and loops back to q0:
//
//	CoeffsToSlots (top) -> EvalMod -> the AES round -> SlotsToCoeffs -> q0 -> ModUp -> top again
func WHTranscipheringParams(logN int, cfg WHChain) (ckks.Parameters, bootstrapping.Parameters, error) {

	// x mod 1 comes first: its depth says how many primes EvalMod eats in the chain.
	mod1Lit := mod1.ParametersLiteral{
		LogScale:        WHLogQ0,
		Mod1Type:        mod1.CosDiscrete,
		Mod1Degree:      cfg.Mod1Degree,
		DoubleAngle:     cfg.DoubleAngle,
		K:               cfg.K,
		LogMessageRatio: WHLogMessageRatio,
	}

	// Written bottom-up, which is the reverse of the order the circuit consumes it.
	logQ := []int{WHLogQ0} // q0, where ScaleDown leaves the state
	for i := 0; i < cfg.S2CDepth; i++ {
		logQ = append(logQ, 60) // SlotsToCoeffs, ct x pt: the prime is free, take a wide one
	}
	roundBottom := len(logQ) - 1                // the level the round has to reach
	logQ = append(logQ, WHLogScale)             // the fused M_r'
	logQ = append(logQ, WHLogScale, WHLogScale) // XOR tree over the four T_r'
	logQ = append(logQ, WHLogScale)             // AddRoundKey
	if cfg.Clean {
		logQ = append(logQ, WHLogScale, WHLogScale) // cleaning
	}
	roundTop := len(logQ) - 1
	for i := 0; i < mod1Lit.Depth(); i++ {
		logQ = append(logQ, WHLogQ0) // EvalMod
	}
	for i := 0; i < cfg.C2SDepth; i++ {
		logQ = append(logQ, WHLogQ0) // CoeffsToSlots
	}

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            logN,
		LogQ:            logQ,
		LogP:            []int{61, 61, 61, 61},
		LogDefaultScale: WHLogScale,
		Xs:              ring.Ternary{H: 256},
	})
	if err != nil {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("params: %w", err)
	}

	// EvalMod reads the ciphertext at scale 2^LogScale to turn m + q0*I into I + m/2^r, so that
	// scale must be q0 rounded to a power of two. Off by one bit and either the reduction lands on
	// 2*(I + eps), which is not what the polynomial approximates, or qDiv drops below 1 and part of
	// the division by the message ratio is pushed into CoeffsToSlots.
	if got := int(math.Round(math.Log2(float64(params.Q()[0])))); mod1Lit.LogScale != got {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf(
			"EvalMod reads at 2^%d but q0 is 2^%d", mod1Lit.LogScale, got)
	}

	// SlotsToCoeffs runs first and ends on q0. Its Scaling is left nil: the evaluator fills it with
	// DefaultScale/(ScalingFactor/MessageRatio), which is exactly the factor ScaleDown then expects.
	S2CParams := dft.MatrixLiteral{
		Type:     dft.HomomorphicDecode,
		LogSlots: params.LogMaxSlots(),
		LevelQ:   cfg.S2CDepth,
		LevelP:   params.MaxLevelP(),
		Levels:   ones(cfg.S2CDepth),
	}

	// CoeffsToSlots runs right after the ModUp, from the top of the chain, and splits the two
	// halves of the coefficient vector into two real-valued ciphertexts.
	C2SParams := dft.MatrixLiteral{
		Type:     dft.HomomorphicEncode,
		Format:   dft.SplitRealAndImag,
		LogSlots: params.LogMaxSlots(),
		LevelQ:   params.MaxLevel(),
		LevelP:   params.MaxLevelP(),
		Levels:   ones(cfg.C2SDepth),
	}

	mod1Lit.LevelQ = params.MaxLevel() - C2SParams.Depth(true)
	if mod1Lit.LevelQ-mod1Lit.Depth() != roundTop {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf(
			"the refresh returns the state at level %d, the round starts at %d",
			mod1Lit.LevelQ-mod1Lit.Depth(), roundTop)
	}
	if roundTop-RoundDepth(cfg) != roundBottom {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf(
			"the round ends at level %d, SlotsToCoeffs starts at %d",
			roundTop-RoundDepth(cfg), roundBottom)
	}

	btpParams := bootstrapping.Parameters{
		ResidualParameters:      params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: S2CParams,
		Mod1ParametersLiteral:   mod1Lit,
		CoeffsToSlotsParameters: C2SParams,
		EphemeralSecretWeight:   32,
		CircuitOrder:            bootstrapping.DecodeThenModUp,
	}
	return params, btpParams, nil
}

// RoundDepth is how many levels one AES middle round consumes: the fused M_r', two for the XOR
// tree, AddRoundKey, and the cleaning when it is on.
func RoundDepth(cfg WHChain) int {
	if cfg.Clean {
		return 6
	}
	return 4
}

// WHRoundTop is the level the refresh returns the state at, i.e. where a round starts.
func WHRoundTop(cfg WHChain) int { return cfg.S2CDepth + RoundDepth(cfg) }

// ones is n levels of one prime each, the Levels form of a DFT matrix literal.
func ones(n int) []int {
	out := make([]int, n)
	for i := range out {
		out[i] = 1
	}
	return out
}
