package cleaning

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Cleaning refines a bit via p(x) = 3x^2 - 2x^3 (p(0+e)=0+e', p(1+e")=1+e*), with e',e* < e, e"
// consuming ceil(log2(deg+1)) = 2 levels.
func Cleaning(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	poly := bignum.NewPolynomial(bignum.Monomial, []complex128{0, 0, 3, -2}, nil)
	polyEval := polynomial.NewEvaluator(params, eval)

	W := ct.Scale
	out, err := polyEval.Evaluate(ct.CopyNew(), poly, W)
	if err != nil {
		return nil, fmt.Errorf("Cleaning: %w", err)
	}
	return out, nil
}
