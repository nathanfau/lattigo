// Package bbbts implements Algorithm 1 (IntRootBoot) of [BKSS24].
// It is a "slim"-order CKKS bootstrap whose EvalMod step is
// replaced by a direct evaluation of the root of unity exp(2*pi*i*m/t): given a
// ciphertext packing integers m_s in {0, ..., t-1} (t = 2^k), IntRootBoot returns a
// ciphertext encrypting the t-th roots of unity {exp(2*pi*i*m_s/t)}, from which the
// bits of each m_s can be recovered (see the bitbatching package).
package bbbts

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/trigo"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const (
	// evalExpR is the number of double-angle squarings of the EvalExp step.
	evalExpR = 3
	// evalExpDeg is the Chebyshev degree of the EvalExp step.
	evalExpDeg = 31
)

// IntRootBoot runs the slim bootstrap pipeline
//
//	SlotsToCoeffs -> ScaleDown -> ModUp -> CoeffsToSlots -> EvalExp
func IntRootBoot(eval *bootstrapping.Evaluator, params ckks.Parameters, ct *rlwe.Ciphertext, k int) (*rlwe.Ciphertext, error) {
	t := 1 << k

	// SlotsToCoeffs (homomorphic decoding); the imaginary part is nil (real input).
	ctSTC, err := eval.SlotsToCoeffs(ct, nil)
	if err != nil {
		return nil, fmt.Errorf("IntRootBoot SlotsToCoeffs: %w", err)
	}

	// ScaleDown to q0/|m|.
	ctSD, _, err := eval.ScaleDown(ctSTC)
	if err != nil {
		return nil, fmt.Errorf("IntRootBoot ScaleDown: %w", err)
	}

	// ModUp (raise the modulus from q0 to qL).
	ctMU, err := eval.ModUp(ctSD)
	if err != nil {
		return nil, fmt.Errorf("IntRootBoot ModUp: %w", err)
	}

	// CoeffsToSlots (homomorphic encoding); only the real part is kept.
	ctReal, _, err := eval.CoeffsToSlots(ctMU)
	if err != nil {
		return nil, fmt.Errorf("IntRootBoot CoeffsToSlots: %w", err)
	}

	// EvalExp replaces EvalMod: reset the scale to the mod-1 scaling factor, while evaluating
	// exp(2*pi*i*x/period) with period = 1/t so the output is exp(2*pi*i*m/t).
	ctReal.Scale = eval.Mod1Parameters.ScalingFactor()
	period := 1.0 / float64(t)

	ctExp, err := trigo.EvalExpDirect(params, eval.Evaluator, ctReal, 1, period, evalExpR, evalExpDeg)
	if err != nil {
		return nil, fmt.Errorf("IntRootBoot EvalExpDirect: %w", err)
	}
	return ctExp, nil
}
