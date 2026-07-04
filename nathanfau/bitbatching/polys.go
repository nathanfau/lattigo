package bitbatching

import (
	"math"
)

type DensePoly []complex128

// ComputePkl computes the coefficients of P_{k,l} by iDFT.
func ComputePkl(k, iter int) DensePoly {
	t := 1 << k
	poly := make(DensePoly, t)
	for n := 0; n < t; n++ {
		var sum complex128
		for m := 0; m < t; m++ {
			if (m>>iter)&1 == 1 {
				angle := -2 * math.Pi * float64(m*n) / float64(t)
				sum += complex(math.Cos(angle), math.Sin(angle))
			}
		}
		poly[n] = sum / complex(float64(t), 0)
	}
	return poly
}

// ComputeQkl extracts Q_{k,l} from P_{k,l}.
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
