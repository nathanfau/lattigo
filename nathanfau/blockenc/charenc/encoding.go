package charenc

import (
	"fmt"
	"math"
	"math/cmplx"
)

type Codec interface {
	Spec() BlockSpec
	EncodeValue(x int) []complex128
	DecodeValue(coords []complex128) (int, error)
	BasisMatrix() [][]complex128
	// Interpolate returns (A, b) with A*EncodeValue(m)+b == values[m] for
	// every m in 0..t-1; each values[m] has length target.Slots.
	Interpolate(values [][]complex128, target BlockSpec) (A [][]complex128, b []complex128, err error)
}

func NewCodec(spec BlockSpec) (Codec, error) {
	switch spec.Alphabet.Kind {
	case BRU:
		return NewBRU(spec.Alphabet.Modulus, spec.Reduced)
	case LBRU:
		return NewLBRU(spec.Alphabet.Modulus, spec.Alphabet.Generator, spec.Reduced)
	case WH:
		return NewWH(spec.Alphabet.Modulus, spec.Reduced)
	case IND:
		return NewIND(spec.Alphabet.Modulus, spec.Reduced, spec.Alphabet.Omitted)
	default:
		return nil, fmt.Errorf("charenc: unsupported encoding kind %v", spec.Alphabet.Kind)
	}
}

// ---------------- BRU ----------------

type bruCodec struct {
	spec BlockSpec
	// start: 1 for reduced (omit trivial character), 0 for non-reduced.
	start int
}

func NewBRU(t int, reduced bool) (Codec, error) {
	if t < 2 {
		return nil, fmt.Errorf("charenc/BRU: modulus must be >= 2, got %d", t)
	}
	slots := t - 1
	start := 1
	if !reduced {
		slots = t
		start = 0
	}
	return &bruCodec{
		spec: BlockSpec{
			Alphabet: Alphabet{Modulus: t, Kind: BRU},
			Reduced:  reduced,
			Slots:    slots,
		},
		start: start,
	}, nil
}

func (c *bruCodec) Spec() BlockSpec { return c.spec }

func (c *bruCodec) EncodeValue(x int) []complex128 {
	t := c.spec.Alphabet.Modulus
	m := positiveMod(x, t)
	out := make([]complex128, c.spec.Slots)
	for j := 0; j < c.spec.Slots; j++ {
		k := c.start + j
		out[j] = rootOfUnity(t, k*m)
	}
	return out
}

func (c *bruCodec) DecodeValue(coords []complex128) (int, error) {
	if len(coords) != c.spec.Slots {
		return 0, fmt.Errorf("charenc/BRU: expected %d coordinates, got %d", c.spec.Slots, len(coords))
	}
	t := c.spec.Alphabet.Modulus
	idx := 1 - c.start
	z := coords[idx]
	theta := math.Atan2(imag(z), real(z))
	m := int(math.Round(theta * float64(t) / (2 * math.Pi)))
	return positiveMod(m, t), nil
}

func (c *bruCodec) BasisMatrix() [][]complex128 { return buildBasisMatrix(c) }

func (c *bruCodec) Interpolate(values [][]complex128, target BlockSpec) ([][]complex128, []complex128, error) {
	return interpolateAffine(c, values, target)
}

// ---------------- LBRU ----------------

type lbruCodec struct {
	spec BlockSpec
	g    int
	dlog []int // dlog[m] = ell with g^ell == m (mod t) for m in 1..t-1
}

// NewLBRU: t must be prime. generator == 0 picks the smallest primitive
// root of Z_t^*. Reduced form: t-1 slots (nonzero indicator then nontrivial
// characters). Non-reduced form prepends a constant slot, giving t slots.
func NewLBRU(t, generator int, reduced bool) (Codec, error) {
	if !isPrime(t) {
		return nil, fmt.Errorf("charenc/LBRU: modulus must be prime, got %d", t)
	}
	if t < 2 {
		return nil, fmt.Errorf("charenc/LBRU: modulus must be >= 2, got %d", t)
	}
	if generator == 0 {
		g, err := smallestPrimitiveRoot(t)
		if err != nil {
			return nil, err
		}
		generator = g
	} else {
		factors := distinctPrimeFactors(t - 1)
		for _, q := range factors {
			if intPowMod(generator, (t-1)/q, t) == 1 {
				return nil, fmt.Errorf("charenc/LBRU: %d is not a primitive root mod %d", generator, t)
			}
		}
	}
	dlog := make([]int, t)
	for i := range dlog {
		dlog[i] = -1
	}
	power := 1 % t
	for ell := 0; ell < t-1; ell++ {
		dlog[power] = ell
		power = (power * generator) % t
	}
	slots := t - 1
	if !reduced {
		slots = t
	}
	return &lbruCodec{
		spec: BlockSpec{
			Alphabet: Alphabet{Modulus: t, Kind: LBRU, Generator: generator},
			Reduced:  reduced,
			Slots:    slots,
		},
		g:    generator,
		dlog: dlog,
	}, nil
}

func (c *lbruCodec) Spec() BlockSpec { return c.spec }

func (c *lbruCodec) EncodeValue(x int) []complex128 {
	t := c.spec.Alphabet.Modulus
	m := positiveMod(x, t)
	out := make([]complex128, c.spec.Slots)
	offset := 0
	if !c.spec.Reduced {
		out[0] = 1
		offset = 1
	}
	if m == 0 {
		return out
	}
	ell := c.dlog[m]
	for k := 0; k < t-1; k++ {
		out[offset+k] = rootOfUnity(t-1, k*ell)
	}
	return out
}

func (c *lbruCodec) DecodeValue(coords []complex128) (int, error) {
	if len(coords) != c.spec.Slots {
		return 0, fmt.Errorf("charenc/LBRU: expected %d coordinates, got %d", c.spec.Slots, len(coords))
	}
	t := c.spec.Alphabet.Modulus
	offset := 0
	if !c.spec.Reduced {
		offset = 1
	}
	if cmplx.Abs(coords[offset]) < 0.5 {
		return 0, nil
	}
	if t == 2 {
		return 1, nil
	}
	z := coords[offset+1]
	theta := math.Atan2(imag(z), real(z))
	ell := int(math.Round(theta * float64(t-1) / (2 * math.Pi)))
	ell = positiveMod(ell, t-1)
	return intPowMod(c.g, ell, t), nil
}

func (c *lbruCodec) BasisMatrix() [][]complex128 { return buildBasisMatrix(c) }

func (c *lbruCodec) Interpolate(values [][]complex128, target BlockSpec) ([][]complex128, []complex128, error) {
	return interpolateAffine(c, values, target)
}

// ---------------- WH ----------------

type whCodec struct {
	spec  BlockSpec
	start int
}

func NewWH(t int, reduced bool) (Codec, error) {
	if _, ok := log2Exact(t); !ok || t < 2 {
		return nil, fmt.Errorf("charenc/WH: modulus must be a power of two >= 2, got %d", t)
	}
	slots := t - 1
	start := 1
	if !reduced {
		slots = t
		start = 0
	}
	return &whCodec{
		spec: BlockSpec{
			Alphabet: Alphabet{Modulus: t, Kind: WH},
			Reduced:  reduced,
			Slots:    slots,
		},
		start: start,
	}, nil
}

func (c *whCodec) Spec() BlockSpec { return c.spec }

func (c *whCodec) EncodeValue(x int) []complex128 {
	t := c.spec.Alphabet.Modulus
	m := positiveMod(x, t)
	out := make([]complex128, c.spec.Slots)
	for j := 0; j < c.spec.Slots; j++ {
		s := c.start + j
		parity := popcount(s&m) & 1
		if parity == 0 {
			out[j] = 1
		} else {
			out[j] = -1
		}
	}
	return out
}

func (c *whCodec) DecodeValue(coords []complex128) (int, error) {
	if len(coords) != c.spec.Slots {
		return 0, fmt.Errorf("charenc/WH: expected %d coordinates, got %d", c.spec.Slots, len(coords))
	}
	t := c.spec.Alphabet.Modulus
	k, _ := log2Exact(t)
	m := 0
	for i := 0; i < k; i++ {
		s := 1 << i
		idx := s - c.start
		bit := 0
		if real(coords[idx]) < 0 {
			bit = 1
		}
		m |= bit << i
	}
	return m, nil
}

func (c *whCodec) BasisMatrix() [][]complex128 { return buildBasisMatrix(c) }

func (c *whCodec) Interpolate(values [][]complex128, target BlockSpec) ([][]complex128, []complex128, error) {
	return interpolateAffine(c, values, target)
}

func popcount(x int) int {
	count := 0
	for x > 0 {
		count += x & 1
		x >>= 1
	}
	return count
}

// ---------------- IND ----------------

type indCodec struct {
	spec BlockSpec
	// values[j] is the alphabet value at slot j; reduced IND skips Omitted.
	values []int
	// slotOf[m] = slot index for value m, or -1 if omitted.
	slotOf []int
}

func NewIND(t int, reduced bool, omitted int) (Codec, error) {
	if t < 1 {
		return nil, fmt.Errorf("charenc/IND: modulus must be >= 1, got %d", t)
	}
	if reduced {
		if omitted < 0 || omitted >= t {
			return nil, fmt.Errorf("charenc/IND: omitted %d out of range [0,%d)", omitted, t)
		}
	}
	slots := t
	if reduced {
		slots = t - 1
	}
	values := make([]int, 0, slots)
	slotOf := make([]int, t)
	for i := range slotOf {
		slotOf[i] = -1
	}
	for m := 0; m < t; m++ {
		if reduced && m == omitted {
			continue
		}
		slotOf[m] = len(values)
		values = append(values, m)
	}
	return &indCodec{
		spec: BlockSpec{
			Alphabet: Alphabet{Modulus: t, Kind: IND, Omitted: omitted},
			Reduced:  reduced,
			Slots:    slots,
		},
		values: values,
		slotOf: slotOf,
	}, nil
}

func (c *indCodec) Spec() BlockSpec { return c.spec }

func (c *indCodec) EncodeValue(x int) []complex128 {
	t := c.spec.Alphabet.Modulus
	m := positiveMod(x, t)
	out := make([]complex128, c.spec.Slots)
	if j := c.slotOf[m]; j >= 0 {
		out[j] = 1
	}
	return out
}

func (c *indCodec) DecodeValue(coords []complex128) (int, error) {
	if len(coords) != c.spec.Slots {
		return 0, fmt.Errorf("charenc/IND: expected %d coordinates, got %d", c.spec.Slots, len(coords))
	}
	best := -1
	bestVal := math.Inf(-1)
	for j, z := range coords {
		r := real(z)
		if r > bestVal {
			bestVal = r
			best = j
		}
	}
	if c.spec.Reduced {
		if bestVal < 0.5 {
			return c.spec.Alphabet.Omitted, nil
		}
	}
	if best < 0 {
		return 0, fmt.Errorf("charenc/IND: empty coordinate vector")
	}
	return c.values[best], nil
}

func (c *indCodec) BasisMatrix() [][]complex128 { return buildBasisMatrix(c) }

func (c *indCodec) Interpolate(values [][]complex128, target BlockSpec) ([][]complex128, []complex128, error) {
	return interpolateAffine(c, values, target)
}
