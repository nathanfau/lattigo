package bitbatching

import (
	"fmt"
	"math"
	"math/bits"
	"math/cmplx"

	ckkspoly "github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	commonpoly "github.com/tuneinsight/lattigo/v6/circuits/common/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// polyCtx evaluates polynomials in X on ONE shared power basis, through lattigo's own
// Paterson-Stockmeyer (EvaluateFromPowerBasis).
//
// lattigo also does the scale bookkeeping, which is the part not to hand-roll: a ciphertext product
// followed by a rescaling lands on \Delta.\Delta/q, not on \Delta, and ckks.Add takes the INTEGER part of the scale
// ratio, so it silently ignores a mismatch of typically 1.0002 and turns it into a relative error on the sum.
type polyCtx struct {
	eval   *ckks.Evaluator
	params ckks.Parameters
	pe     *ckkspoly.Evaluator
	pb     commonpoly.PowerBasis
	W      rlwe.Scale
	lvIn   int
}

func newPolyCtx(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext) *polyCtx {
	return &polyCtx{
		eval:   eval,
		params: params,
		pe:     ckkspoly.NewEvaluator(params, eval),
		pb:     commonpoly.NewPowerBasis(ct, bignum.Monomial),
		W:      ct.Scale,
		lvIn:   ct.Level(),
	}
}

// evalAt evaluates coeffs (index i for X^i) and returns a ciphertext carrying exactly target.
// A nil result means the polynomial is zero.
func (c *polyCtx) evalAt(coeffs []complex128, target rlwe.Scale) (*rlwe.Ciphertext, error) {
	cc := chop(coeffs)
	if polyDegree(cc) < 0 {
		return nil, nil
	}
	poly := bignum.NewPolynomial(bignum.Monomial, cc, nil)
	out, err := c.pe.EvaluateFromPowerBasis(c.pb, poly, target)
	if err != nil {
		return nil, fmt.Errorf("EvaluateFromPowerBasis: %w", err)
	}
	if err = c.settle(out, target); err != nil {
		return nil, err
	}
	return out, nil
}

// power returns X^n from the shared basis, generating it if needed.
func (c *polyCtx) power(n int) (*rlwe.Ciphertext, error) {
	if err := c.pb.GenPower(n, false, c.eval); err != nil {
		return nil, fmt.Errorf("GenPower X^%d: %w", n, err)
	}
	return c.pb.Value[n], nil
}

// settle closes the gap when the last rescaling lags: this fork's polynomial evaluator can return a
// prime above the target.
func (c *polyCtx) settle(ct *rlwe.Ciphertext, target rlwe.Scale) error {
	for i := 0; i < 2 && ct.Scale.Div(target).Float64() > 1.5; i++ {
		if err := c.eval.Rescale(ct, ct); err != nil {
			return fmt.Errorf("settle Rescale: %w", err)
		}
	}
	return checkScale("polynomial result", ct.Scale, target)
}

// normDepth adds to a polynomial's own depth the level settle spends closing the lagging rescale.
// That level is not always taken, so pinLevel forces it: it is what makes normTarget's arithmetic
// exact rather than accidentally right.
func normDepth(depth int) int { return depth + 1 }

// normTarget is the scale a polynomial must be produced at so that its product with u, which
// happens at level mulLvl, lands exactly on target: (T . s_u) / q[mulLvl] = target.
func (c *polyCtx) normTarget(depth int, u *rlwe.Ciphertext, target rlwe.Scale) rlwe.Scale {
	mulLvl := min(c.lvIn-normDepth(depth), u.Level())
	return target.Mul(qScale(c.params, mulLvl)).Div(u.Scale)
}

// pinLevel drops ct to exactly depth levels below the input, the level normTarget assumed for the
// later product. A ciphertext sitting LOWER than predicted is left alone: checkScale will say so.
func (c *polyCtx) pinLevel(ct *rlwe.Ciphertext, depth int) {
	if d := ct.Level() - (c.lvIn - normDepth(depth)); d > 0 {
		c.eval.DropLevel(ct, d)
	}
}

// polyDepth is the number of levels the evaluator spends on a polynomial of this degree. It chops
// first, like the evaluation does: the inverse-DFT dust would report degree t-1 for every bit
// function instead of the real support.
func polyDepth(coeffs []complex128) int {
	d := polyDegree(chop(coeffs))
	if d < 1 {
		return 0
	}
	return bits.Len(uint(d))
}

// polyDegree is the highest index carrying a non-zero coefficient, -1 for the zero polynomial.
// Not bignum.Polynomial.Degree, which is len(Coeffs)-1 and so reports t-1 on these sparse arrays.
func polyDegree(coeffs []complex128) int {
	for i := len(coeffs) - 1; i >= 0; i-- {
		if coeffs[i] != 0 {
			return i
		}
	}
	return -1
}

// chop zeroes the coefficients that are only inverse-DFT rounding dust. Without it a 1e-17 in an
// even-only support defeats lattigo's parity detection and doubles the basis it builds.
func chop(p []complex128) []complex128 {
	maxAbs := 0.0
	for _, c := range p {
		maxAbs = math.Max(maxAbs, cmplx.Abs(c))
	}
	tol := 1e-9 * math.Max(1, maxAbs)
	out := make([]complex128, len(p))
	for i, c := range p {
		if cmplx.Abs(c) > tol {
			out[i] = c
		}
	}
	return out
}

func qScale(params ckks.Parameters, lvl int) rlwe.Scale {
	return utils.RescalingFactor(params, lvl)
}

// checkScale fails loudly when a ciphertext does not carry the scale the caller is about to assume
// it carries. A silent divergence would show up as a precision floor and nothing else.
func checkScale(what string, got, want rlwe.Scale) error {
	if !got.InDelta(want, 40) {
		return fmt.Errorf("%s: scale 2^%.2f off the target", what, -got.Log2Delta(want))
	}
	return nil
}
