package cleaning

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// The monomial coefficients of the polynomials below, lowest degree first.
var (
	basicCoeffs        = []complex128{0, 0, 3, -2}
	smootherCoeffs     = []complex128{0, 0, 0, 10, -15, 6}
	verySmootherCoeffs = []complex128{0, 0, 0, 0, 35, -84, 70, -20}
	signCoeffs         = []complex128{0, 1.5, 0, -0.5}
)

// evalAt evaluates the polynomial on ct and lands the result on the scale W, exactly: lattigo's
// polynomial evaluator picks its constants for that target, whatever the input scale.
func evalAt(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext, coeffs []complex128,
	W rlwe.Scale, name string) (*rlwe.Ciphertext, error) {
	poly := bignum.NewPolynomial(bignum.Monomial, coeffs, nil)
	out, err := polynomial.NewEvaluator(params, eval).Evaluate(ct.CopyNew(), poly, W)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", name, err)
	}
	return out, nil
}

// Cleaning refines a bit via p(x) = 3x^2 - 2x^3 (p(0+e)=0+e', p(1+e")=1+e*), with e',e* < e, e"
// consuming ceil(log2(deg+1)) = 2 levels.
func Cleaning(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return evalAt(params, eval, ct, basicCoeffs, ct.Scale, "Cleaning")
}

// SmootherCleaning refines a bit via p(x) = 6x^5-15x^4+10x^3
// consuming ceil(log2(deg+1)) = 3 levels.
func SmootherCleaning(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return evalAt(params, eval, ct, smootherCoeffs, ct.Scale, "Cleaning")
}

// VerySmootherCleaning refines a bit via p(x) = -20x^7+70x^6-84x^5+35x^4
// consuming ceil(log2(deg+1)) = 3 levels.
func VerySmootherCleaning(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return evalAt(params, eval, ct, verySmootherCoeffs, ct.Scale, "Cleaning")
}

// SignCleaning refines a slot carrying +1 or -1.
// It consumes the same ceil(log2(deg+1)) = 2 levels as Cleaning.
func SignCleaning(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return evalAt(params, eval, ct, signCoeffs, ct.Scale, "SignCleaning")
}
