// Package trigo evaluates trigonometric functions (cos, sin, exp) homomorphically
// via Chebyshev approximation (Paterson-Stockmeyer) on the interval [-K, K].
//
// The base frequency is fT = period * 2^r: the polynomials approximate the function
// at frequency fT, and r "double-angle" squarings bring the result to the target
// frequency 1/period. EvalCos and EvalSin leave the squarings to the caller (used by
// Algo1 in the real / conjugate-invariant context); EvalExpDirect applies them itself
// (used by the standard / complex context).
package trigo

import (
	"fmt"
	"math"
	"math/big"
	"math/cmplx"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// EvalCos evaluates cos(2*pi*x / fT) with fT = period * 2^r on [-K, K], via a
// degree-'degree' real Chebyshev approximation. The r squarings are NOT applied here.
// Cost: 1 (change of basis, skipped if K == 1) + ceil(log2(degree+1)) (poly).
func EvalCos(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext, K int, period float64, r, degree int) (*rlwe.Ciphertext, error) {
	fT := period * math.Pow(2, float64(r))
	poly := chebyPolyReal(float64(K), degree, func(x float64) float64 {
		return math.Cos(2 * math.Pi * x / fT)
	})
	return evalCheby(params, eval, ct, poly, "EvalCos")
}

// EvalSin is the sinus-variant of EvalCos
func EvalSin(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext, K int, period float64, r, degree int) (*rlwe.Ciphertext, error) {
	fT := period * math.Pow(2, float64(r))
	poly := chebyPolyReal(float64(K), degree, func(x float64) float64 {
		return math.Sin(2 * math.Pi * x / fT)
	})
	return evalCheby(params, eval, ct, poly, "EvalSin")
}

// EvalExpDirect evaluates exp(2*pi*i*x / period) on [-K, K], via a degree-'degree'
// complex Chebyshev approximation at the base frequency fT = period * 2^r, then applies
// r squarings (double-angle) to reach the target frequency 1/period.
// Cost: 1 (change of basis, skipped if K == 1) + ceil(log2(degree+1)) (poly) + r (squarings).
func EvalExpDirect(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext, K int, period float64, r, degree int) (*rlwe.Ciphertext, error) {
	fT := period * math.Pow(2, float64(r))
	poly := chebyPolyComplex(float64(K), degree, func(x complex128) complex128 {
		return cmplx.Exp(complex(0, 2*math.Pi*real(x)/fT))
	})
	out, err := evalCheby(params, eval, ct, poly, "EvalExpDirect")
	if err != nil {
		return nil, err
	}
	for i := 0; i < r; i++ {
		if err := eval.MulRelin(out, out, out); err != nil {
			return nil, fmt.Errorf("EvalExpDirect squaring %d: %w", i+1, err)
		}
		if err := eval.Rescale(out, out); err != nil {
			return nil, fmt.Errorf("EvalExpDirect Rescale squaring %d: %w", i+1, err)
		}
	}
	return out, nil
}

// Exp is evaluated directly rather than as cos + i*sin: benchmarks (not kept here) showed
// the same precision either way, but the direct route is about 80% faster (one complex-polynomial
// evaluation instead of two real-polynomial).

// evalCheby applies the Chebyshev change of basis (ct' = scalar*ct + constant) then
// evaluates poly on the result. The change of basis is skipped when it is the identity
// (scalar == 1 and constant == 0, e.g. the interval [-1, 1]).
func evalCheby(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext, poly bignum.Polynomial, name string) (*rlwe.Ciphertext, error) {
	scalar, constant := poly.ChangeOfBasis()

	work := ct.CopyNew()
	if math.Abs(bigToF64(scalar)-1) > 1e-9 || math.Abs(bigToF64(constant)) > 1e-9 {
		if err := eval.Mul(work, scalar, work); err != nil {
			return nil, fmt.Errorf("%s Mul CoB: %w", name, err)
		}
		if err := eval.Add(work, constant, work); err != nil {
			return nil, fmt.Errorf("%s Add CoB: %w", name, err)
		}
		if err := eval.Rescale(work, work); err != nil {
			return nil, fmt.Errorf("%s Rescale CoB: %w", name, err)
		}
	}

	polyEval := polynomial.NewEvaluator(params, eval)
	out, err := polyEval.Evaluate(work, poly, work.Scale)
	if err != nil {
		return nil, fmt.Errorf("%s Evaluate: %w", name, err)
	}
	return out, nil
}

// bigToF64 converts a *big.Float to a float64 (losing precision).
func bigToF64(x *big.Float) float64 {
	f, _ := x.Float64()
	return f
}

// chebyPolyReal builds the degree-'degree' real Chebyshev approximation of f on [-K, K].
func chebyPolyReal(K float64, degree int, f func(float64) float64) bignum.Polynomial {
	const prec uint = 128
	return bignum.ChebyshevApproximation(f, bignum.Interval{
		A:     *bignum.NewFloat(-K, prec),
		B:     *bignum.NewFloat(K, prec),
		Nodes: degree,
	})
}

// chebyPolyComplex builds the degree-'degree' complex Chebyshev approximation of f on [-K, K].
func chebyPolyComplex(K float64, degree int, f func(complex128) complex128) bignum.Polynomial {
	const prec uint = 128
	return bignum.ChebyshevApproximation(f, bignum.Interval{
		A:     *bignum.NewFloat(-K, prec),
		B:     *bignum.NewFloat(K, prec),
		Nodes: degree,
	})
}
