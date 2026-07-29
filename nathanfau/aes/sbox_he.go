// All the homomorphic S-box (SubByte) variants and their helpers are in this file.
// It is purely an organizational split from aes_he.go
//
// The three variants (V2, V3, V4) share one implementation, subByte, and differ ONLY by how
// aggressively they defer the relin + rescale of the built ANF monomials.
package aes

import (
	"fmt"
	"math/bits"
	"sync"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
)

// sboxANF is the cleartext ANF description of the AES S-box: it depends on the S-box only, so it is
// computed once and shared by every SubByte call.
type sboxANF struct {
	coeffs  [8][256]int  // coeffs[i][S] = multilinear coefficient a_S of (sbox>>i)&1
	needed  []int        // monomial masks S != 0 with a non-zero coefficient in some output bit
	factors map[int]bool // monomials reused as a half in the balanced build of another monomial
}

// anf returns the (once-computed) ANF description of the AES S-box.
var anf = sync.OnceValue(func() *sboxANF {
	coeffs := sboxCoeffs()
	needed := neededMonomials(coeffs)
	return &sboxANF{coeffs: coeffs, needed: needed, factors: factorMonomials(needed)}
})

// sboxCoeffs[i][S] = multilinear coefficient a_S of (sbox>>i)&1.
func sboxCoeffs() [8][256]int {
	var c [8][256]int
	for i := 0; i < 8; i++ {
		var a [256]int
		for T := 0; T < 256; T++ {
			a[T] = int((sbox[T] >> uint(i)) & 1)
		}
		for j := 0; j < 8; j++ {
			for S := 0; S < 256; S++ {
				if S>>uint(j)&1 == 1 {
					a[S] -= a[S^(1<<uint(j))]
				}
			}
		}
		c[i] = a
	}
	return c
}

// neededMonomials returns the ANF monomials (subset masks S != 0) that appear with a non-zero
// coefficient in at least one output bit of the S-box.
func neededMonomials(coeffs [8][256]int) []int {
	var needed []int
	var seen [256]bool
	for i := 0; i < 8; i++ {
		for S := 1; S < 256; S++ {
			if coeffs[i][S] != 0 && !seen[S] {
				seen[S] = true
				needed = append(needed, S)
			}
		}
	}
	return needed
}

// balancedSplit splits a monomial mask S into two halves of cardinality floor(k/2) and ceil(k/2)
// (e.g. x0x1x2x3 = (x0x1)(x2x3), not (x0x1x2)(x3)), giving multiplicative depth 3 for k <= 8.
func balancedSplit(S int) (loMask, hiMask int) {
	var idx []int
	for j := 0; j < 8; j++ {
		if S>>uint(j)&1 == 1 {
			idx = append(idx, j)
		}
	}
	h := len(idx) / 2
	for _, j := range idx[:h] {
		loMask |= 1 << uint(j)
	}
	for _, j := range idx[h:] {
		hiMask |= 1 << uint(j)
	}
	return loMask, hiMask
}

// factorMonomials returns the set of monomials (masks, degree >= 2) that are REUSED as a half in
// the balanced build of some needed monomial. A factor is consumed by a further multiplication, so
// it can never be kept lazy: it must be relinearized + rescaled at build time.
func factorMonomials(needed []int) map[int]bool {
	factors := map[int]bool{}
	seen := map[int]bool{}
	for j := 0; j < 8; j++ {
		seen[1<<uint(j)] = true
	}
	var walk func(S int)
	walk = func(S int) {
		if seen[S] {
			return
		}
		seen[S] = true
		lo, hi := balancedSplit(S)
		if bits.OnesCount(uint(lo)) >= 2 {
			factors[lo] = true
		}
		if bits.OnesCount(uint(hi)) >= 2 {
			factors[hi] = true
		}
		walk(lo)
		walk(hi)
	}
	for _, S := range needed {
		walk(S)
	}
	return factors
}

// lazyFn reports whether a monomial mask is kept LAZY, i.e. left as a degree-2 product at scale
// ~Delta^2 (utils.MulLeveledLazy: no relin, no rescale) instead of being reduced back to degree 1
// at build time (utils.MulLeveled).
type lazyFn func(S int) bool

// lazyAbove builds the laziness predicate of a variant: a monomial is lazy when it is a LEAF (not
// reused as a build factor, so nothing multiplies it again) and its degree is >= minDeg.
// minDeg = 0 disables laziness entirely.
func (s *sboxANF) lazyAbove(minDeg int) lazyFn {
	return func(S int) bool {
		return minDeg > 0 && !s.factors[S] && bits.OnesCount(uint(S)) >= minDeg
	}
}

// buildMonomials builds every needed monomial from the 8 input bits by balanced splitting, taking
// the lazy path (degree-2, deferred relin + rescale) for the monomials isLazy selects and the
// standard leveled path (one relin + one rescale each) for the others.
func (a *Evaluator) buildMonomials(inp ByteHE, s *sboxANF, isLazy lazyFn) (map[int]*rlwe.Ciphertext, error) {
	mono := make(map[int]*rlwe.Ciphertext, len(s.needed)+8)
	for j := 0; j < 8; j++ {
		mono[1<<uint(j)] = inp[j].CopyNew()
	}
	var build func(S int) (*rlwe.Ciphertext, error)
	build = func(S int) (*rlwe.Ciphertext, error) {
		if m, ok := mono[S]; ok {
			return m, nil
		}
		lo, hi := balancedSplit(S)
		l, err := build(lo)
		if err != nil {
			return nil, err
		}
		h, err := build(hi)
		if err != nil {
			return nil, err
		}
		var m *rlwe.Ciphertext
		if isLazy(S) {
			m, err = utils.MulLeveledLazy(a.eval, l, h) // defer relin + rescale (degree 2)
		} else {
			m, err = utils.MulLeveled(a.eval, l, h) // relin + rescale now (degree 1)
		}
		if err != nil {
			return nil, fmt.Errorf("buildMonomials S=%d: %w", S, err)
		}
		mono[S] = m
		return m, nil
	}
	for _, S := range s.needed {
		if _, err := build(S); err != nil {
			return nil, err
		}
	}
	return mono, nil
}

// sboxSum computes each output bit as the ANF-weighted sum of the monomials, which consumes no
// level: only Add/Sub and integer Mul by constants. All monomials are first brought to a common
// level (Add/Sub match the residual scale drift), then the terms are split into two accumulators:
// the eager monomials (degree-1, scale ~Delta) and the lazy ones (degree-2, scale ~Delta^2). The
// lazy sum is relinearized + rescaled ONCE per output bit before the two are combined, which is
// what makes the lazy variants cheap.
func (a *Evaluator) sboxSum(mono map[int]*rlwe.Ciphertext, s *sboxANF, isLazy lazyFn) (out ByteHE, err error) {
	monoCts := make([]*rlwe.Ciphertext, 0, len(s.needed))
	for _, S := range s.needed {
		monoCts = append(monoCts, mono[S])
	}
	utils.FlattenLevels(a.eval, monoCts)

	addTerm := func(accPtr **rlwe.Ciphertext, S, c int) error {
		acc := *accPtr
		switch {
		case c == 1:
			if acc == nil {
				acc = mono[S].CopyNew()
			} else if e := a.eval.Add(acc, mono[S], acc); e != nil {
				return e
			}
		case c == -1:
			if acc == nil {
				acc = mono[S].CopyNew()
				if e := a.eval.Mul(acc, -1, acc); e != nil {
					return e
				}
			} else if e := a.eval.Sub(acc, mono[S], acc); e != nil {
				return e
			}
		default:
			t := mono[S].CopyNew()
			if e := a.eval.Mul(t, c, t); e != nil {
				return e
			}
			if acc == nil {
				acc = t
			} else if e := a.eval.Add(acc, t, acc); e != nil {
				return e
			}
		}
		*accPtr = acc
		return nil
	}

	for i := 0; i < 8; i++ {
		var accE, accL *rlwe.Ciphertext
		for _, S := range s.needed {
			c := s.coeffs[i][S]
			if c == 0 {
				continue
			}
			acc, kind := &accE, "eager"
			if isLazy(S) {
				acc, kind = &accL, "lazy"
			}
			if err = addTerm(acc, S, c); err != nil {
				return out, fmt.Errorf("sboxSum bit %d %s S=%d: %w", i, kind, S, err)
			}
		}
		if accL != nil {
			if err = a.eval.Relinearize(accL, accL); err != nil {
				return out, fmt.Errorf("sboxSum bit %d relin: %w", i, err)
			}
			utils.Ops.Relin++
			if err = a.eval.Rescale(accL, accL); err != nil {
				return out, fmt.Errorf("sboxSum bit %d rescale: %w", i, err)
			}
			utils.Ops.Rescale++
		}
		var acc *rlwe.Ciphertext
		switch {
		case accE != nil && accL != nil:
			utils.AlignLevels(a.eval, accE, accL)
			if err = a.eval.Add(accE, accL, accE); err != nil {
				return out, fmt.Errorf("sboxSum bit %d combine: %w", i, err)
			}
			acc = accE
		case accE != nil:
			acc = accE
		case accL != nil:
			acc = accL
		default:
			return out, fmt.Errorf("sboxSum bit %d: empty sum", i)
		}
		if a0 := s.coeffs[i][0]; a0 != 0 {
			if err = a.eval.Add(acc, a0, acc); err != nil {
				return out, fmt.Errorf("sboxSum bit %d add const: %w", i, err)
			}
		}
		out[i] = acc
	}
	return out, nil
}

// subByte applies the AES S-box to a bit-sliced encrypted byte (8 ciphertexts) through its ANF: it
// builds the needed monomials from the input bits, then sums them per output bit.
func (a *Evaluator) subByte(inp ByteHE, lazyMinDeg int) (ByteHE, error) {
	s := anf()
	isLazy := s.lazyAbove(lazyMinDeg)
	mono, err := a.buildMonomials(inp, s, isLazy)
	if err != nil {
		return ByteHE{}, err
	}
	return a.sboxSum(mono, s, isLazy)
}

// SubByteV2 is the non-lazy [BCKK25] baseline: every monomial is relinearized + rescaled at build
// time. Cost: one relin + one rescale per monomial of degree >= 2, i.e. 247 = 255 - 8 per byte.
func (a *Evaluator) SubByteV2(inp ByteHE) (ByteHE, error) { return a.subByte(inp, 0) }

// SubByteV3 is fully lazy: only the 61 monomials reused as a build factor are relinearized +
// rescaled, the 186 leaf monomials stay degree-2 at scale ~Delta^2 and their relin + rescale is
// deferred to one pair per output-bit accumulator. Cost: 61 + 8 = 69 relin and 69 rescale per byte.
func (a *Evaluator) SubByteV3(inp ByteHE) (ByteHE, error) { return a.subByte(inp, 2) }

// SubByteV4 is a less aggressive lazy strategy than SubByteV3: only the leaf monomials of degree
// >= 4 are kept lazy, the factors and the degree-2/3 leaves are relinearized + rescaled at build
// time. Cost: 61 factors + 29 low-degree leaves + 8 accumulators = 98 relin and 98 rescale per byte.
func (a *Evaluator) SubByteV4(inp ByteHE) (ByteHE, error) { return a.subByte(inp, 4) }

func (a *Evaluator) SubByte(inp ByteHE, version int) (ByteHE, error) {
	switch version {
	case 2:
		return a.SubByteV2(inp)
	case 3:
		return a.SubByteV3(inp)
	case 4:
		return a.SubByteV4(inp)
	default:
		return ByteHE{}, fmt.Errorf("SubByte: unknown version %d (want 2..4)", version)
	}
}
