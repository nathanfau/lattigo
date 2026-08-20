package whaes

import "fmt"

// The round matrices are Kronecker products of a 16x16 pattern over the AES bytes and a 256x256
// block acting inside one byte, built separately so a single A_alpha serves every fusion.

// Identity is I_n.
func Identity(n int) [][]complex128 {
	m := make([][]complex128, n)
	for i := range m {
		m[i] = make([]complex128, n)
		m[i][i] = 1
	}
	return m
}

// Zeros is the r x c zero matrix.
func Zeros(r, c int) [][]complex128 {
	m := make([][]complex128, r)
	for i := range m {
		m[i] = make([]complex128, c)
	}
	return m
}

// Kron is the Kronecker product, (a tensor b)[i*m+k][j*n+l] = a[i][j] * b[k][l]: the byte pattern
// indexes the blocks, the inner matrix the slots inside a block, which is the WH packing.
func Kron(a, b [][]complex128) [][]complex128 {
	p := len(a)
	if p == 0 {
		return nil
	}
	q := len(a[0])
	m := len(b)
	if m == 0 {
		return nil
	}
	n := len(b[0])

	out := Zeros(p*m, q*n)
	for i := 0; i < p; i++ {
		for j := 0; j < q; j++ {
			v := a[i][j]
			if v == 0 {
				continue
			}
			for k := 0; k < m; k++ {
				rowOut := out[i*m+k]
				rowB := b[k]
				for l := 0; l < n; l++ {
					rowOut[j*n+l] = v * rowB[l]
				}
			}
		}
	}
	return out
}

// Add sums matrices of the same shape, how M_r' assembles from its three alpha terms.
func Add(ms ...[][]complex128) ([][]complex128, error) {
	if len(ms) == 0 {
		return nil, fmt.Errorf("whaes.Add: no matrix")
	}
	r, c := len(ms[0]), len(ms[0][0])
	out := Zeros(r, c)
	for x, m := range ms {
		if len(m) != r || len(m[0]) != c {
			return nil, fmt.Errorf("whaes.Add: matrix %d is %dx%d, want %dx%d", x, len(m), len(m[0]), r, c)
		}
		for i := 0; i < r; i++ {
			for j := 0; j < c; j++ {
				out[i][j] += m[i][j]
			}
		}
	}
	return out, nil
}

// Mul composes the byte patterns, P_{r',alpha}.P_pi, before they are tensored with A_alpha.
func Mul(a, b [][]complex128) ([][]complex128, error) {
	if len(a) == 0 || len(b) == 0 {
		return nil, fmt.Errorf("whaes.Mul: empty matrix")
	}
	inner := len(a[0])
	if inner != len(b) {
		return nil, fmt.Errorf("whaes.Mul: %d columns against %d rows", inner, len(b))
	}
	out := Zeros(len(a), len(b[0]))
	for i := range a {
		for k := 0; k < inner; k++ {
			v := a[i][k]
			if v == 0 {
				continue
			}
			for j := range b[k] {
				out[i][j] += v * b[k][j]
			}
		}
	}
	return out, nil
}
