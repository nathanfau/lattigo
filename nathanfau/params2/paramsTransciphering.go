// Package params holds the CKKS and bootstrapping parameter sets for the transciphering tests
// (nathanfau/transciphering): the k=4 Algo1 pipeline (TranscipheringParams). The other packages'
// tests build their own parameters.
package params2

import (
	"fmt"
	"math"
	"math/big"
	"slices"
	"strconv"
	"strings"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
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

// TranscipheringParams is TranscipheringParamsWith at the default cleaning depth, 2 levels, i.e.
// cleaning.Cleaning, the default extraction, bitbatching.BitExtract at k levels, and the pipeline's
// own shape. It must stay in step with transciphering.DefaultCleanDepth, which cannot be imported
// here: transciphering already depends on this package.
func TranscipheringParams(logN, k int) (ckks.Parameters, bootstrapping.Parameters, error) {
	return TranscipheringParamsWith(logN, k, 2, k, Shape{})
}

// DefaultLogP is the auxiliary modulus: 8 primes of 42 bits, i.e. dnum = Ceil(#Q/8) = 4 digits on
// the 31-prime chain. 42 is the smallest size still covering the widest digit (326 bits), and
// anything above that margin is log(QP) spent on noise already below the rounding.
func DefaultLogP() []int { return []int{42, 42, 42, 42, 42, 42, 42, 42} }

// DefaultLogSTC is what SlotsToCoeffs spends: one 60-bit prime, hence one level. It is what makes
// the bottom key-switching digit 22 bits wider than a run of 38-bit primes; splitting it over more
// primes spreads the same budget over more levels, and lifts the whole circuit by as many.
func DefaultLogSTC() []int { return []int{60} }

// STCLevels is the SlotsToCoeffs depth a given set of primes buys: one level each.
func STCLevels(logSTC []int) []int {
	lv := make([]int, len(logSTC))
	for i := range lv {
		lv[i] = 1
	}
	return lv
}

// LogScale is the default scale of the pipeline, and the size of every prime of the chain but q0,
// the SlotsToCoeffs block and the zone (Shape.LogZone).
const LogScale = 38

// DefaultLogQ0 is the bottom prime. Lattigo couples it to the scale by q0 = LogScale + k for a
// message ratio of 2^k; we drop the coupling and run at LogScale, ratio 1, the 1/t being carried by
// the SlotsToCoeffs scaling instead (see the README) -- measured, the ratio costs nothing.
const DefaultLogQ0 = LogScale

// Shape is the part of the parameter set the sweeps vary. The zero value is the pipeline's own.
type Shape struct {
	LogP   []int // sizes of the P primes; nil = DefaultLogP
	LogSTC []int // sizes of the SlotsToCoeffs primes, one level each; nil = DefaultLogSTC
	LogQ0  int   // size of the bottom prime; 0 = the chain scale plus DefaultLogQ0 - LogScale, the
	// default message ratio. The ratio follows it: LogMessageRatio = LogQ0 - LogQi, so with 38-bit
	// primes 42 gives 2^4 and 38 gives 1.

	// LogQi is the size of the chain's primes, every Q prime but q0 and the SlotsToCoeffs ones, and
	// the scale the pipeline runs at: 0 = LogScale. Neither P nor the SlotsToCoeffs primes follow
	// it, so moving it moves the widest key-switching digit against a fixed P.
	LogQi int

	// LogZone is the size of the primes from Conv_{Real->Cplx} up to the bit extraction, i.e. the
	// AES circuit and the bottom of the refresh; 0 = LogQi, no zone. Those stages multiply
	// ciphertexts together, so the scale they run at must be the size of their primes: LogZone is
	// both, and ZoneScale is that scale. The bootstrap above and below keeps LogQi.
	LogZone int
}

// logQi is LogQi with its zero resolved.
func (sh Shape) logQi() int {
	if sh.LogQi == 0 {
		return LogScale
	}
	return sh.LogQi
}

// logZone is LogZone with its zero resolved: the chain's own primes, no zone.
func (sh Shape) logZone() int {
	if sh.LogZone == 0 {
		return sh.logQi()
	}
	return sh.LogZone
}

// HasZone reports whether sh runs a zone, i.e. primes of another size than the chain's.
func (sh Shape) HasZone() bool { return sh.logZone() != sh.logQi() }

// ZoneScale is the scale a ciphertext must carry inside the zone sh describes: 2^LogZone, which
// is the pipeline's own default scale when there is no zone.
func (sh Shape) ZoneScale() rlwe.Scale { return rlwe.NewScale(math.Exp2(float64(sh.logZone()))) }

// ParseLogSTC reads a SlotsToCoeffs shape the way the test flags write it: the prime sizes, "30,30",
// or n equal primes, "2x30".
func ParseLogSTC(s string) ([]int, error) {
	s = strings.TrimSpace(s)
	if n, b, ok := strings.Cut(s, "x"); ok {
		count, err1 := strconv.Atoi(strings.TrimSpace(n))
		bits, err2 := strconv.Atoi(strings.TrimSpace(b))
		if err1 != nil || err2 != nil || count < 1 || bits < 1 {
			return nil, fmt.Errorf("SlotsToCoeffs shape %q: want <count>x<bits>, e.g. 2x30", s)
		}
		return slices.Repeat([]int{bits}, count), nil
	}
	var logSTC []int
	for _, f := range strings.Split(s, ",") {
		b, err := strconv.Atoi(strings.TrimSpace(f))
		if err != nil || b < 1 {
			return nil, fmt.Errorf("SlotsToCoeffs shape %q: want prime sizes, e.g. 30,30, or 2x30", s)
		}
		logSTC = append(logSTC, b)
	}
	return logSTC, nil
}

// FormatLogSTC writes a shape the way ParseLogSTC reads it, "2x30" when its primes are all equal.
func FormatLogSTC(logSTC []int) string {
	if len(logSTC) > 1 && slices.Min(logSTC) == slices.Max(logSTC) {
		return fmt.Sprintf("%dx%d", len(logSTC), logSTC[0])
	}
	s := make([]string, len(logSTC))
	for i, b := range logSTC {
		s[i] = strconv.Itoa(b)
	}
	return strings.Join(s, ",")
}

// TranscipheringParamsWith is the Algo1 pipeline parameter set (nathanfau/transciphering) on a
// fully specified shape, sh's zero fields meaning the pipeline's own. cleanDepth is the cleaning
// polynomial's depth and extractLv the bit extraction's, k or k+1; both lengthen the chain above
// the circuit, where a longer SlotsToCoeffs lengthens it below and shifts every level up.
func TranscipheringParamsWith(logN, k, cleanDepth, extractLv int, sh Shape) (ckks.Parameters, bootstrapping.Parameters, error) {
	logP, logSTC, logQ0, logQi, zone := sh.LogP, sh.LogSTC, sh.LogQ0, sh.LogQi, sh.logZone()
	if logP == nil {
		logP = DefaultLogP()
	}
	if logSTC == nil {
		logSTC = DefaultLogSTC()
	}
	if logQi == 0 {
		logQi = LogScale
	}
	if logQi < 1 {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("logQi %d, want a prime size in bits", logQi)
	}
	if zone < 1 || zone > logQi {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("logZone %d, want between 1 and the chain's %d", zone, logQi)
	}
	if logQ0 == 0 {
		logQ0 = logQi + DefaultLogQ0 - LogScale
	}
	if logQ0 < logQi || logQ0 > logQi+k {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("logQ0 %d, want between the chain scale %d (message ratio 1) and %d (ratio 2^k)", logQ0, logQi, logQi+k)
	}
	if len(logSTC) == 0 {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("logSTC is empty: SlotsToCoeffs spends at least one prime")
	}
	if len(logP) == 0 {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("logP is empty: the pipeline key-switches, so it needs at least one P prime")
	}
	if cleanDepth < 2 || cleanDepth > 3 {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("cleanDepth %d, want 2 or 3", cleanDepth)
	}
	if extractLv < k || extractLv > k+1 {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("extractLv %d, want %d or %d", extractLv, k, k+1)
	}
	logQ := []int{logQ0}
	// qi repeats the chain prime n times, zq the zone's.
	qi := func(n int) []int { return slices.Repeat([]int{logQi}, n) }
	zq := func(n int) []int { return slices.Repeat([]int{zone}, n) }
	logQ = append(logQ, logSTC...) // SlotsToCoeffs, one level per prime
	// The zone: every stage from here to the bit extraction runs at ZoneScale.
	logQ = append(logQ, zq(1)...)          // Conv_{Real->Cplx}
	logQ = append(logQ, zq(3)...)          // SubBytes
	logQ = append(logQ, zq(cleanDepth)...) // cleaning
	logQ = append(logQ, zq(1)...)          // AddRoundKey
	logQ = append(logQ, zq(3)...)          // MixColumns
	// refresh
	logQ = append(logQ, zq(1)...)         // Conv_{Cplx->Real}
	logQ = append(logQ, zq(extractLv)...) // bit extraction
	// End of the zone.
	logQ = append(logQ, qi(3)...) // 3 levels for squaring
	logQ = append(logQ, qi(1)...) // extractExp
	logQ = append(logQ, qi(1)...) // Conv_{Real->Cplx}
	logQ = append(logQ, qi(5)...) // EvalCos
	logQ = append(logQ, qi(1)...) // Conv_{Cplx->Real}
	logQ = append(logQ, qi(3)...) // CoeffsToSlots

	params, err := ciThenStd(logN, logQ, logP, logQi)
	if err != nil {
		return ckks.Parameters{}, bootstrapping.Parameters{}, err
	}

	S2CParams := dft.MatrixLiteral{
		Type:     dft.HomomorphicDecode,
		LogSlots: params.LogMaxSlots(),
		LevelP:   params.MaxLevelP(),
		Levels:   STCLevels(logSTC),
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
		LogScale:        logQi,
		Mod1Type:        mod1.CosDiscrete,
		Mod1Degree:      2 * ((1 << k) - 1),
		K:               1 << k,
		LogMessageRatio: logQ0 - logQi,
	}
	mod1P, err := mod1.NewParametersFromLiteral(params, Mod1Params)
	if err != nil {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("mod1 params: %w", err)
	}
	// SlotsToCoeffs hands the mod-1 circuit a message already divided by t, whatever q0 is. The
	// message ratio carries that division only when q0 = LogScale + k, so the 1/t is written here:
	// the same number in that case, and the right one in both.
	F := mod1P.ScalingFactor().Float64() / (params.DefaultScale().Float64() * float64(uint(1)<<k))
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
