package bitbatching

// DensePoly is a polynomial by its coefficients, index i holding the one of X^i.
type DensePoly []complex128

// ComputePkl returns the interpolant P_{k,l} of the l-th bit over Z_{2^k}, i.e. the polynomial of
// degree < 2^k that maps omega^m to (m >> l) & 1 on the 2^k-th roots of unity. By Lemma 3 it is
// supported on 0 and on the odd multiples of 2^{k-l-1}.
func ComputePkl(k, l int) DensePoly {
	return Lagrange(BitTargets(k, l))
}

// ComputeQkl extracts the LOWER half of the spectrum of P_{k,l}, as the coefficients of a dense
// polynomial in Y = X^{2^{k-l}}: q[p] is the coefficient of X^{step.(2p+1)}, step = 2^{k-l-1}.
//
// Half is enough because the bit is real, so the upper half of the spectrum is the conjugate of the
// lower one and comes back from a single conjugation. Returns nil for l = 0, whose support is
// {0, t/2} and which has no lower half at all.
func ComputeQkl(k, iter int, pkl DensePoly) DensePoly {
	if iter == 0 {
		return nil
	}
	numTerms := 1 << (iter - 1)
	q := make(DensePoly, numTerms)
	step := 1 << (k - iter - 1)
	for p := 0; p < numTerms; p++ {
		q[p] = pkl[step*(2*p+1)]
	}
	return q
}
