package charenc

import (
	"fmt"
	"math"
	"math/cmplx"
	"sync"
)

const DefaultDecodeTolerance = 1e-6

func positiveMod(x, m int) int {
	r := x % m
	if r < 0 {
		r += m
	}
	return r
}

func rootOfUnity(t, k int) complex128 {
	if t <= 0 {
		return complex(math.NaN(), math.NaN())
	}
	theta := 2 * math.Pi * float64(positiveMod(k, t)) / float64(t)
	return complex(math.Cos(theta), math.Sin(theta))
}

func intPowMod(base, exp, modulus int) int {
	if modulus == 1 {
		return 0
	}
	result := 1 % modulus
	b := positiveMod(base, modulus)
	for exp > 0 {
		if exp&1 == 1 {
			result = (result * b) % modulus
		}
		exp >>= 1
		b = (b * b) % modulus
	}
	return result
}

func isPrime(n int) bool {
	if n < 2 {
		return false
	}
	if n%2 == 0 {
		return n == 2
	}
	for d := 3; d*d <= n; d += 2 {
		if n%d == 0 {
			return false
		}
	}
	return true
}

func isPowerOfTwo(n int) bool {
	return n > 0 && (n&(n-1)) == 0
}

func log2Exact(n int) (int, bool) {
	if !isPowerOfTwo(n) {
		return 0, false
	}
	k := 0
	for n > 1 {
		n >>= 1
		k++
	}
	return k, true
}

func distinctPrimeFactors(n int) []int {
	out := []int{}
	if n <= 1 {
		return out
	}
	x := n
	for d := 2; d*d <= x; d++ {
		if x%d == 0 {
			out = append(out, d)
			for x%d == 0 {
				x /= d
			}
		}
	}
	if x > 1 {
		out = append(out, x)
	}
	return out
}

// smallestPrimitiveRoot requires p prime, p >= 2.
func smallestPrimitiveRoot(p int) (int, error) {
	if !isPrime(p) {
		return 0, fmt.Errorf("smallestPrimitiveRoot: %d is not prime", p)
	}
	if p == 2 {
		return 1, nil
	}
	factors := distinctPrimeFactors(p - 1)
	for g := 2; g < p; g++ {
		ok := true
		for _, q := range factors {
			if intPowMod(g, (p-1)/q, p) == 1 {
				ok = false
				break
			}
		}
		if ok {
			return g, nil
		}
	}
	return 0, fmt.Errorf("smallestPrimitiveRoot: no generator found for %d", p)
}

func buildBasisMatrix(c Codec) [][]complex128 {
	spec := c.Spec()
	t := spec.Alphabet.Modulus
	rows := make([][]complex128, t)
	for m := 0; m < t; m++ {
		rows[m] = c.EncodeValue(m)
	}
	return rows
}

// SolveSquareSystem solves M*x = y. For repeated solves with the same M,
// reuse LUFactorize(M).Solve.
func SolveSquareSystem(M [][]complex128, y []complex128) ([]complex128, error) {
	lu, err := LUFactorize(M)
	if err != nil {
		return nil, err
	}
	return lu.Solve(y)
}

// LUFactorization packs L (unit diagonal implicit) and U into one row-major
// matrix; perm is the partial-pivoting row order.
type LUFactorization struct {
	lu   [][]complex128
	perm []int
	n    int
}

func LUFactorize(M [][]complex128) (*LUFactorization, error) {
	n := len(M)
	if n == 0 {
		return nil, fmt.Errorf("LUFactorize: empty matrix")
	}
	lu := make([][]complex128, n)
	for i := range lu {
		if len(M[i]) != n {
			return nil, fmt.Errorf("LUFactorize: row %d has %d cols, expected %d", i, len(M[i]), n)
		}
		row := make([]complex128, n)
		copy(row, M[i])
		lu[i] = row
	}
	perm := make([]int, n)
	for i := range perm {
		perm[i] = i
	}
	for col := 0; col < n; col++ {
		piv := col
		bestAbs := cmplx.Abs(lu[col][col])
		for r := col + 1; r < n; r++ {
			if a := cmplx.Abs(lu[r][col]); a > bestAbs {
				bestAbs = a
				piv = r
			}
		}
		if bestAbs < 1e-12 {
			return nil, fmt.Errorf("LUFactorize: singular at column %d", col)
		}
		if piv != col {
			lu[col], lu[piv] = lu[piv], lu[col]
			perm[col], perm[piv] = perm[piv], perm[col]
		}
		invPiv := 1 / lu[col][col]
		for r := col + 1; r < n; r++ {
			factor := lu[r][col] * invPiv
			lu[r][col] = factor
			if factor == 0 {
				continue
			}
			luCol := lu[col]
			luR := lu[r]
			for j := col + 1; j < n; j++ {
				luR[j] -= factor * luCol[j]
			}
		}
	}
	return &LUFactorization{lu: lu, perm: perm, n: n}, nil
}

func (f *LUFactorization) Solve(y []complex128) ([]complex128, error) {
	if len(y) != f.n {
		return nil, fmt.Errorf("LUFactorization.Solve: y length %d, expected %d", len(y), f.n)
	}
	b := make([]complex128, f.n)
	for i, p := range f.perm {
		b[i] = y[p]
	}
	for i := 1; i < f.n; i++ {
		luI := f.lu[i]
		s := b[i]
		for j := 0; j < i; j++ {
			s -= luI[j] * b[j]
		}
		b[i] = s
	}
	x := make([]complex128, f.n)
	for i := f.n - 1; i >= 0; i-- {
		luI := f.lu[i]
		s := b[i]
		for j := i + 1; j < f.n; j++ {
			s -= luI[j] * x[j]
		}
		x[i] = s / luI[i]
	}
	return x, nil
}

var luCache sync.Map

// LUFactorizeCached computes the factorization once per key; the caller must
// choose a key that uniquely identifies M.
func LUFactorizeCached(key string, M [][]complex128) (*LUFactorization, error) {
	if v, ok := luCache.Load(key); ok {
		return v.(*LUFactorization), nil
	}
	lu, err := LUFactorize(M)
	if err != nil {
		return nil, err
	}
	if v, loaded := luCache.LoadOrStore(key, lu); loaded {
		return v.(*LUFactorization), nil
	}
	return lu, nil
}

// interpolateAffine returns (A, b) with A*EncodeValue(m)+b == values[m].
// Reduced encodings solve the augmented basis [1|M]; the first coefficient
// becomes the bias. Non-reduced encodings carry the constant coordinate in M
// already, so b is zero.
func interpolateAffine(c Codec, values [][]complex128, target BlockSpec) ([][]complex128, []complex128, error) {
	spec := c.Spec()
	t := spec.Alphabet.Modulus
	if len(values) != t {
		return nil, nil, fmt.Errorf("interpolateAffine: got %d input rows, expected %d", len(values), t)
	}
	outSlots := target.Slots
	for i, row := range values {
		if len(row) != outSlots {
			return nil, nil, fmt.Errorf("interpolateAffine: row %d has %d cols, target has %d", i, len(row), outSlots)
		}
	}

	M := c.BasisMatrix()
	var sysMat [][]complex128
	withBias := spec.Reduced
	if withBias {
		sysMat = make([][]complex128, t)
		for i := 0; i < t; i++ {
			row := make([]complex128, 1+spec.Slots)
			row[0] = 1
			copy(row[1:], M[i])
			sysMat[i] = row
		}
	} else {
		sysMat = M
	}
	if len(sysMat) > 0 && len(sysMat[0]) != t {
		return nil, nil, fmt.Errorf("interpolateAffine: basis is not square (%d x %d) for alphabet size %d", len(sysMat), len(sysMat[0]), t)
	}

	A := make([][]complex128, outSlots)
	for j := range A {
		A[j] = make([]complex128, spec.Slots)
	}
	b := make([]complex128, outSlots)

	cacheKey := fmt.Sprintf("interpolateAffine/%s/m=%d/r=%v/g=%d/o=%d/withBias=%v",
		spec.Alphabet.Kind, spec.Alphabet.Modulus, spec.Reduced,
		spec.Alphabet.Generator, spec.Alphabet.Omitted, withBias)
	lu, err := LUFactorizeCached(cacheKey, sysMat)
	if err != nil {
		return nil, nil, fmt.Errorf("interpolateAffine: factorize: %w", err)
	}

	for j := 0; j < outSlots; j++ {
		yj := make([]complex128, t)
		for i := 0; i < t; i++ {
			yj[i] = values[i][j]
		}
		coef, err := lu.Solve(yj)
		if err != nil {
			return nil, nil, fmt.Errorf("interpolateAffine: solve j=%d: %w", j, err)
		}
		if withBias {
			b[j] = coef[0]
			copy(A[j], coef[1:])
		} else {
			copy(A[j], coef)
		}
	}
	return A, b, nil
}
