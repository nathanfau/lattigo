package bitbatching

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// The three ways of getting the k bits of m back out of a ciphertext whose slots hold omega^m,
// omega = exp(2i.pi/2^k). Same signature, same output order, same output scale and level, so a
// caller can hold any of them behind one variable -- transciphering.Config.Extract does.
//
//	BitExtract        half spectrum, [BBBTS] 3.2       k levels    error LINEAR in the input error
//	BitExtractInterp  plain interpolation, Lemma 3     k levels    error linear, worse constant
//	BitExtractClean   interpolation + cleaning, 3.3    k+1 levels  error QUADRATIC
//
// BitExtractInterp is the control, called by no one else: same interpolants as BitExtract but
// whole, same targets as BitExtractClean but without the bivariate lift, so a comparison against it
// attributes a gain to the right cause.
//
// Below them, the evaluation of h ([BBBTS] 3.3) that BitExtractClean drives. Only the real-target
// form h = F0 + 2.Re(U + x.conj(x).W) has an evaluator here, which is all the bit functions need;
// its coefficients are built in cleaninterp.go.

// BitExtract homomorphically extracts the k bits of m from ct encrypting omega^m, omega =
// exp(2i.pi/2^k).
//
// [BBBTS] 3.2. P_{k,l} is supported on the odd multiples of step = 2^{k-l-1} (Lemma 3), and the bit
// being real, its upper spectrum is the conjugate of its lower one. Only HALF is evaluated and one
// conjugation brings the rest back, constant term included. The error stays LINEAR: (t/4).eps at
// worst, on the LSB.
//
// The LSB is a special case: X^{t/2} is (-1)^m, already real on the roots, so it needs neither a
// half spectrum nor a conjugation.
func BitExtract(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext, k int) ([]*rlwe.Ciphertext, error) {
	return BitExtractTo(params, eval, ct, k, ct.Scale)
}

// BitExtractTo is BitExtract with the bits landed on the scale W instead of the input's. W must be
// close to ct.Scale; the landing is exact and costs no level.
func BitExtractTo(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext, k int, W rlwe.Scale) ([]*rlwe.Ciphertext, error) {
	if k < 1 {
		return nil, fmt.Errorf("BitExtract: k must be >= 1")
	}

	t := 1 << k
	c := newPolyCtx(params, eval, ct, W)

	results := make([]*rlwe.Ciphertext, k)
	for l := 0; l < k; l++ {
		pkl := ComputePkl(k, l)
		step := 1 << (k - l - 1)

		coeffs := make([]complex128, t)
		if l == 0 {
			coeffs[0], coeffs[t/2] = pkl[0], pkl[t/2]
		} else {
			for p, v := range ComputeQkl(k, l, pkl) {
				coeffs[step*(2*p+1)] = v
			}
		}

		out, err := c.evalAt(coeffs, c.W)
		if err != nil {
			return nil, fmt.Errorf("BitExtract bit %d: %w", l, err)
		}
		if out == nil {
			return nil, fmt.Errorf("BitExtract bit %d: empty interpolant", l)
		}

		if l > 0 {
			if out, err = realFold(eval, pkl[0], out); err != nil {
				return nil, fmt.Errorf("BitExtract bit %d fold: %w", l, err)
			}
		}
		results[l] = out
	}

	utils.FlattenLevels(eval, results)
	return results, nil
}

// BitExtractInterp extracts the k bits through the PLAIN Lagrange interpolation of Lemma 3, no
// cleaning: the control. Error LINEAR, with a slightly worse constant than the half spectrum.
//
// Same k levels as the other two: its degree t-1 and BitExtract's X^{t/2} round up to the same
// ceil(log2(.)). The half spectrum buys ciphertext products, not depth.
func BitExtractInterp(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext, k int) ([]*rlwe.Ciphertext, error) {
	if k < 1 {
		return nil, fmt.Errorf("BitExtractInterp: k must be >= 1")
	}

	c := newPolyCtx(params, eval, ct, ct.Scale)

	results := make([]*rlwe.Ciphertext, k)
	for l := 0; l < k; l++ {
		out, err := c.evalAt(ComputePkl(k, l), c.W)
		if err != nil {
			return nil, fmt.Errorf("BitExtractInterp bit %d: %w", l, err)
		}
		if out == nil {
			return nil, fmt.Errorf("BitExtractInterp bit %d: empty interpolant", l)
		}
		results[l] = out
	}

	utils.FlattenLevels(eval, results)
	return results, nil
}

// BitExtractClean extracts the k bits like BitExtract, but through the combined interpolation and
// cleaning, so the bits come out cleaned. Every target is real, so the symmetrised rank t/2
// collapses h to two polynomials per bit instead of four.
//
// k+1 levels, one more than BitExtract (5 at k = 4, measured): U costs k on its own, the x.conj(x)
// side is a level shallower but pays it back on the product with u, and one more on the rescaling
// settle closes -- that last one is the only one that could still be argued with.
func BitExtractClean(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext, k int) ([]*rlwe.Ciphertext, error) {
	return BitExtractCleanTo(params, eval, ct, k, ct.Scale)
}

// BitExtractCleanTo is BitExtractClean with the bits landed on the scale W instead of the input's,
// as BitExtractTo.
func BitExtractCleanTo(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext, k int, W rlwe.Scale) ([]*rlwe.Ciphertext, error) {
	if k < 1 {
		return nil, fmt.Errorf("BitExtractClean: k must be >= 1")
	}

	c := newPolyCtx(params, eval, ct, W)

	u, err := rootNorm(eval, ct, c.W)
	if err != nil {
		return nil, fmt.Errorf("BitExtractClean: %w", err)
	}

	results := make([]*rlwe.Ciphertext, k)
	for l := 0; l < k; l++ {
		ci, err := NewCleanInterp(BitTargets(k, l), true)
		if err != nil {
			return nil, fmt.Errorf("BitExtractClean bit %d: %w", l, err)
		}
		if !ci.Symmetric() {
			return nil, fmt.Errorf("BitExtractClean bit %d: V != U, the real-target collapse does not hold", l)
		}

		var out *rlwe.Ciphertext
		if l == 0 {
			out, err = cleanLSB(c, ci, u)
		} else {
			out, err = evalSym(c, ci, u)
		}
		if err != nil {
			return nil, fmt.Errorf("BitExtractClean bit %d: %w", l, err)
		}
		if err = checkScale(fmt.Sprintf("BitExtractClean bit %d result", l), out.Scale, c.W); err != nil {
			return nil, err
		}
		out.Scale = c.W
		results[l] = out
	}

	utils.FlattenLevels(eval, results)
	return results, nil
}

// evalSym evaluates h in its real-target form, F0 + 2.Re(U + u.W), the body BitExtractClean runs on
// every bit but the LSB. U comes out on W, and W is produced at the scale that makes u.W land back
// on W, so everything meets there and the additions are exact.
func evalSym(c *polyCtx, ci CleanInterp, u *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	depthW := polyDepth(ci.W)

	ctU, err := c.evalAt(ci.U, c.W)
	if err != nil {
		return nil, fmt.Errorf("U: %w", err)
	}
	ctW, err := c.evalAt(ci.W, c.normTarget(depthW, u, c.W))
	if err != nil {
		return nil, fmt.Errorf("W: %w", err)
	}
	c.pinLevel(ctW, depthW)

	return assembleSym(c.eval, ci.F0, ctU, ctW, u, c.W)
}

// cleanLSB is h for the LSB, whose U and W both reduce to the monomial X^{t/2}. Evaluated as two
// polynomials, u.W would sit a level below U and set the depth of the whole extraction; factoring
// as h = F0 + 2.Re((A + C.u).X^{t/2}) keeps it level with the other bits. lambda puts the result on
// W, which the polynomial evaluator would otherwise have done.
func cleanLSB(c *polyCtx, ci CleanInterp, u *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	s := ci.T / 2
	xs, err := c.power(s)
	if err != nil {
		return nil, err
	}

	mulLvl := min(u.Level()-1, xs.Level())
	lam := complex(c.W.Mul(qScale(c.params, mulLvl)).Div(u.Scale.Mul(xs.Scale)).Float64(), 0)

	acc := u.CopyNew()
	if err = c.eval.Mul(acc, ci.W[s]*lam, acc); err != nil {
		return nil, fmt.Errorf("cleanLSB Mul C: %w", err)
	}
	if err = c.eval.RescaleTo(acc, c.W, acc); err != nil {
		return nil, fmt.Errorf("cleanLSB RescaleTo C: %w", err)
	}
	if err = c.eval.Add(acc, ci.U[s]*lam, acc); err != nil {
		return nil, fmt.Errorf("cleanLSB Add A: %w", err)
	}

	xc := xs.CopyNew()
	utils.AlignLevels(c.eval, acc, xc)
	if err = c.eval.MulRelin(acc, xc, acc); err != nil {
		return nil, fmt.Errorf("cleanLSB MulRelin: %w", err)
	}
	if err = c.eval.RescaleTo(acc, c.W, acc); err != nil {
		return nil, fmt.Errorf("cleanLSB RescaleTo: %w", err)
	}
	// The product lands on s_u.s_x/q, not on W: lambda corrected the VALUE, the only freedom a
	// scalar multiplication gives. Once that is checked, assigning W is a true statement.
	if err = checkScale("cleanLSB", acc.Scale, u.Scale.Mul(xs.Scale).Div(qScale(c.params, mulLvl))); err != nil {
		return nil, err
	}
	acc.Scale = c.W

	return realFold(c.eval, ci.F0, acc)
}

// rootNorm returns u = x.conj(x), the |x|^2 factor every trinomial carries on its third monomial.
// It is REAL, which is what lets assembleSym conjugate U + u.W in one go: conj(u.W) = u.conj(W).
func rootNorm(eval *ckks.Evaluator, ct *rlwe.Ciphertext, W rlwe.Scale) (*rlwe.Ciphertext, error) {
	cj := ct.CopyNew()
	if err := eval.Conjugate(cj, cj); err != nil {
		return nil, fmt.Errorf("rootNorm Conjugate: %w", err)
	}
	out := ct.CopyNew()
	if err := eval.MulRelin(out, cj, out); err != nil {
		return nil, fmt.Errorf("rootNorm MulRelin: %w", err)
	}
	if err := eval.RescaleTo(out, W, out); err != nil {
		return nil, fmt.Errorf("rootNorm RescaleTo: %w", err)
	}
	return out, nil
}

// assembleSym builds h = F0 + T + conj(T) with T = U + u.W: one conjugation for the whole thing,
// since u is real.
func assembleSym(eval *ckks.Evaluator, f0 complex128, ctU, ctW, u *rlwe.Ciphertext, W rlwe.Scale) (*rlwe.Ciphertext, error) {
	t, err := mulByNorm(eval, ctW, u, W)
	if err != nil {
		return nil, err
	}
	if t == nil {
		t = ctU
	} else if ctU != nil {
		if err = checkScale("assembleSym U vs u.W", t.Scale, ctU.Scale); err != nil {
			return nil, err
		}
		t.Scale = ctU.Scale
		if err = eval.Add(t, ctU, t); err != nil {
			return nil, fmt.Errorf("assembleSym U + u.W: %w", err)
		}
	}
	if t == nil {
		return nil, fmt.Errorf("assembleSym: h has no non-constant term")
	}
	return realFold(eval, f0, t)
}

// realFold returns f0 + 2.Re(ct). Conjugation is a key switch, moving neither scale nor level, so
// the addition is exact.
func realFold(eval *ckks.Evaluator, f0 complex128, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	out := ct.CopyNew()
	cj := ct.CopyNew()
	if err := eval.Conjugate(cj, cj); err != nil {
		return nil, fmt.Errorf("realFold Conjugate: %w", err)
	}
	if err := eval.Add(out, cj, out); err != nil {
		return nil, fmt.Errorf("realFold Add conj: %w", err)
	}
	if err := eval.Add(out, f0, out); err != nil {
		return nil, fmt.Errorf("realFold Add f0: %w", err)
	}
	return out, nil
}

// mulByNorm returns u.ct, or nil when ct is nil.
func mulByNorm(eval *ckks.Evaluator, ct, u *rlwe.Ciphertext, W rlwe.Scale) (*rlwe.Ciphertext, error) {
	if ct == nil {
		return nil, nil
	}
	out := ct.CopyNew()
	un := u.CopyNew()
	utils.AlignLevels(eval, out, un)
	if err := eval.MulRelin(out, un, out); err != nil {
		return nil, fmt.Errorf("mulByNorm MulRelin: %w", err)
	}
	if err := eval.RescaleTo(out, W, out); err != nil {
		return nil, fmt.Errorf("mulByNorm RescaleTo: %w", err)
	}
	return out, nil
}
