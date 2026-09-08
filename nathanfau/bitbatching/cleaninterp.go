package bitbatching

import (
	"fmt"
	"math"
	"math/cmplx"
)

// Combined interpolation and cleaning for roots of unity ([BBBTS] 3.3), the coefficient side: the
// bivariate h below agrees with the Lagrange interpolant f on the roots, but its terms linear in
// eps cancel, so the error at zeta^j + eps is quadratic -- no h1 step, no depth beyond log t.
// Evaluated by [BitExtractClean]; s = floor(t/2).
//
//	h(x) = f_0
//	     + sum_{k=1..s}   f_k/t ( k.conj(x)^{t-k} + (t-k)(k+1).x^k - k(t-k).x^{k+1}.conj(x) )
//	     + sum_{k=s+1..t-1} f_k/t ( (t-k).x^k + k(t-k+1).conj(x)^{t-k} - k(t-k).conj(x)^{t-k+1}.x )

// Lagrange returns the size-t inverse DFT of y, i.e. the coefficients of the interpolant of degree
// < t through it. Hand-rolled because lattigo has no interpolation helper.
func Lagrange(y []complex128) DensePoly {
	t := len(y)
	f := make(DensePoly, t)
	for k := 0; k < t; k++ {
		var sum complex128
		for j := 0; j < t; j++ {
			angle := -2 * math.Pi * float64(j*k) / float64(t)
			sum += y[j] * complex(math.Cos(angle), math.Sin(angle))
		}
		f[k] = sum / complex(float64(t), 0)
	}
	return f
}

// CleanInterp is h regrouped by monomial, as four polynomials in x indexed by the power of x they
// multiply (index 0 unused, the constant term being F0):
//
//	h(x) = F0 + U(x) + conj(V(x)) + x.conj(x) . ( W(x) + conj(Z(x)) )
type CleanInterp struct {
	T          int        // t, the number of roots of unity
	F0         complex128 // f_0
	U, V, W, Z DensePoly  // length t each
}

// NewCleanInterp builds h for the targets y (y[j] for zeta^j). Two departures from the article:
//
//   - the conjugate side switches at t-s-1, not s. The two differ at the single rank t/2, where
//     taking s leaves h EXACT ON THE ROOTS but kills the cleaning -- error back to order 1, i.e.
//     15.7 bits at t = 16 on the LSB, the one bit function with a t/2 coefficient.
//   - sym averages the two admissible forms of the rank t/2 trinomial, which restores V == U and
//     Z == W on real targets: two polynomial evaluations instead of four. A no-op when f_{t/2} = 0.
func NewCleanInterp(y []complex128, sym bool) (CleanInterp, error) {
	t := len(y)
	if t < 2 {
		return CleanInterp{}, fmt.Errorf("NewCleanInterp: t = %d, want >= 2", t)
	}
	if sym && t%2 != 0 {
		return CleanInterp{}, fmt.Errorf("NewCleanInterp: sym needs t even, got t = %d", t)
	}

	f := Lagrange(y)
	s := t / 2
	thr := t - s - 1
	ct := complex(float64(t), 0)

	// A_k, C_k: coefficients of x^k, plain and in factor of x.conj(x). B_j, D_j: the same on the
	// conj(x) side, contributed by rank k = t-j.
	A := make(DensePoly, t)
	B := make(DensePoly, t)
	C := make(DensePoly, t)
	D := make(DensePoly, t)
	for k := 1; k < t; k++ {
		if k <= s {
			A[k] = f[k] / ct * complex(float64((t-k)*(k+1)), 0)
		} else {
			A[k] = f[k] / ct * complex(float64(t-k), 0)
		}
	}
	for j := 1; j < t; j++ {
		if j <= thr {
			B[j] = f[t-j] / ct * complex(float64((t-j)*(j+1)), 0)
		} else {
			B[j] = f[t-j] / ct * complex(float64(t-j), 0)
		}
	}
	for k := 1; k <= s; k++ {
		C[k] = -f[k] / ct * complex(float64(k*(t-k)), 0)
	}
	for j := 1; j <= thr; j++ {
		D[j] = -f[t-j] / ct * complex(float64((t-j)*j), 0)
	}
	if sym {
		fs := f[s] / complex(float64(2*t), 0)
		A[s] = fs * complex(float64(s*(s+2)), 0)
		B[s] = A[s]
		C[s] = -fs * complex(float64(s*s), 0)
		D[s] = C[s]
	}

	ci := CleanInterp{T: t, F0: f[0], U: A, W: C}
	ci.V = make(DensePoly, t)
	ci.Z = make(DensePoly, t)
	for j := range B {
		ci.V[j] = cmplx.Conj(B[j])
		ci.Z[j] = cmplx.Conj(D[j])
	}
	return ci, nil
}

// Symmetric reports whether V == U and Z == W, so that h collapses to two polynomials.
func (ci CleanInterp) Symmetric() bool {
	tol := 1e-10 * math.Max(1, ci.maxAbs())
	for j := range ci.U {
		if cmplx.Abs(ci.U[j]-ci.V[j]) > tol || cmplx.Abs(ci.W[j]-ci.Z[j]) > tol {
			return false
		}
	}
	return true
}

func (ci CleanInterp) maxAbs() float64 {
	m := 0.0
	for _, p := range []DensePoly{ci.U, ci.V, ci.W, ci.Z} {
		for _, c := range p {
			if a := cmplx.Abs(c); a > m {
				m = a
			}
		}
	}
	return m
}

// BitTargets returns the targets of the l-th bit over Z_{2^k}: y[m] = (m >> l) & 1.
func BitTargets(k, l int) []complex128 {
	t := 1 << k
	y := make([]complex128, t)
	for m := 0; m < t; m++ {
		y[m] = complex(float64((m>>l)&1), 0)
	}
	return y
}
