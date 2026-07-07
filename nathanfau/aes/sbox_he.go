// All the homomorphic S-box (SubByte) variants and their helpers are in this file.
// It is piurely a organizational split from aes_he.go
package aes

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
)

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

// buildMonomials builds every needed monomial from the 8 input bits, a monomial of cardinality k splits into halves of
// cardinality ceil(k/2) and floor(k/2) (e.g. x0x1x2x3 = (x0x1)(x2x3), not (x0x1x2)(x3)), giving
// multiplicative depth 3 for k <= 8.
func (a *Evaluator) buildMonomials(inp ByteHE, needed []int, mul func(x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error)) (map[int]*rlwe.Ciphertext, error) {
	mono := make(map[int]*rlwe.Ciphertext, len(needed)+8)
	for j := 0; j < 8; j++ {
		mono[1<<uint(j)] = inp[j].CopyNew()
	}
	var build func(S int) (*rlwe.Ciphertext, error)
	build = func(S int) (*rlwe.Ciphertext, error) {
		if m, ok := mono[S]; ok {
			return m, nil
		}
		var bits []int
		for j := 0; j < 8; j++ {
			if S>>uint(j)&1 == 1 {
				bits = append(bits, j)
			}
		}
		h := len(bits) / 2
		var loMask, hiMask int
		for _, j := range bits[:h] {
			loMask |= 1 << uint(j)
		}
		for _, j := range bits[h:] {
			hiMask |= 1 << uint(j)
		}
		lo, err := build(loMask)
		if err != nil {
			return nil, err
		}
		hi, err := build(hiMask)
		if err != nil {
			return nil, err
		}
		m, err := mul(lo, hi)
		if err != nil {
			return nil, fmt.Errorf("buildMonomials S=%d: %w", S, err)
		}
		mono[S] = m
		return m, nil
	}
	for _, S := range needed {
		if _, err := build(S); err != nil {
			return nil, err
		}
	}
	return mono, nil
}

// sboxSum computes each output bit as the ANF-weighted sum of the monomials. All monomials must
// already sit at the same level (Add/Sub match the residual scale drift). It consumes no level
// and no relin/rescale: only Add/Sub and integer Mul/Add by constants.
func (a *Evaluator) sboxSum(mono map[int]*rlwe.Ciphertext, needed []int, coeffs [8][256]int) (out ByteHE, err error) {
	for i := 0; i < 8; i++ {
		var acc *rlwe.Ciphertext
		for _, S := range needed {
			c := coeffs[i][S]
			switch {
			case c == 0:
				continue
			case c == 1:
				if acc == nil {
					acc = mono[S].CopyNew()
				} else if err = a.eval.Add(acc, mono[S], acc); err != nil {
					return out, fmt.Errorf("sboxSum bit %d add S=%d: %w", i, S, err)
				}
			case c == -1:
				if acc == nil {
					acc = mono[S].CopyNew()
					if err = a.eval.Mul(acc, -1, acc); err != nil {
						return out, fmt.Errorf("sboxSum bit %d neg S=%d: %w", i, S, err)
					}
				} else if err = a.eval.Sub(acc, mono[S], acc); err != nil {
					return out, fmt.Errorf("sboxSum bit %d sub S=%d: %w", i, S, err)
				}
			default:
				t := mono[S].CopyNew()
				if err = a.eval.Mul(t, c, t); err != nil { // integer cmult: scale unchanged
					return out, fmt.Errorf("sboxSum bit %d cmult S=%d: %w", i, S, err)
				}
				if acc == nil {
					acc = t
				} else if err = a.eval.Add(acc, t, acc); err != nil {
					return out, fmt.Errorf("sboxSum bit %d add-cmult S=%d: %w", i, S, err)
				}
			}
		}
		if acc == nil {
			return out, fmt.Errorf("sboxSum bit %d: empty sum", i)
		}
		if a0 := coeffs[i][0]; a0 != 0 { // constant term a_empty
			if err = a.eval.Add(acc, a0, acc); err != nil {
				return out, fmt.Errorf("sboxSum bit %d add const: %w", i, err)
			}
		}
		out[i] = acc
	}
	return out, nil
}

// SubByteV1 applies the AES S-box to a bit-sliced encrypted byte (8 ciphertexts). It builds
// the ANF monomials with the SQUARING product (utils.MulBits) and flattens levels by squaring
// (utils.FlattenBitLevels), so every monomial stays exactly on the canonical scale chain.
func (a *Evaluator) SubByteV1(inp ByteHE) (ByteHE, error) {
	coeffs := sboxCoeffs()
	needed := neededMonomials(coeffs)
	mono, err := a.buildMonomials(inp, needed, func(x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
		return utils.MulBits(a.eval, x, y)
	})
	if err != nil {
		return ByteHE{}, err
	}
	monoCts := make([]*rlwe.Ciphertext, 0, len(needed))
	for _, S := range needed {
		monoCts = append(monoCts, mono[S])
	}
	if err := utils.FlattenBitLevels(a.eval, monoCts); err != nil {
		return ByteHE{}, fmt.Errorf("SubByte flatten: %w", err)
	}
	return a.sboxSum(mono, needed, coeffs)
}

// SubByteV2 applies the AES S-box like SubByteV1 but builds the monomials the standard LEVELED way
// (utils.MulLeveled) and flattens by DropLevel (utils.FlattenLevels), never squaring for level
// alignment. This removes every alignment squaring: the cost is exactly one relin + one rescale
// per built monomial of degree >= 2 (247 = 255 - 8 for the AES S-box), the non-lazy [BCKK25]
// baseline.
func (a *Evaluator) SubByteV2(inp ByteHE) (ByteHE, error) {
	coeffs := sboxCoeffs()
	needed := neededMonomials(coeffs)
	mono, err := a.buildMonomials(inp, needed, func(x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
		return utils.MulLeveled(a.eval, x, y)
	})
	if err != nil {
		return ByteHE{}, err
	}
	monoCts := make([]*rlwe.Ciphertext, 0, len(needed))
	for _, S := range needed {
		monoCts = append(monoCts, mono[S])
	}
	utils.FlattenLevels(a.eval, monoCts)
	return a.sboxSum(mono, needed, coeffs)
}

// balancedSplit splits a monomial mask S into two halves of cardinality floor(k/2) and
// ceil(k/2)
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
// the balanced build of some needed monomial. In the lazy S-box these are the only monomials
// that must be relinearized+rescaled.
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

// buildMonomialsLazy builds every needed monomial like buildMonomials, but relinearizes+rescales
// only the FACTOR monomials (utils.MulLeveled, degree-1) and keeps the LEAF monomials as
// degree-2 products at scale ~Delta^2 (utils.MulLeveledLazy, no relin/rescale).
func (a *Evaluator) buildMonomialsLazy(inp ByteHE, needed []int, factors map[int]bool) (map[int]*rlwe.Ciphertext, error) {
	mono := make(map[int]*rlwe.Ciphertext, len(needed)+8)
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
		if factors[S] {
			m, err = utils.MulLeveled(a.eval, l, h) // reused: relin + rescale now (degree 1)
		} else {
			m, err = utils.MulLeveledLazy(a.eval, l, h) // leaf: defer relin + rescale (degree 2)
		}
		if err != nil {
			return nil, fmt.Errorf("buildMonomialsLazy S=%d: %w", S, err)
		}
		mono[S] = m
		return m, nil
	}
	for _, S := range needed {
		if _, err := build(S); err != nil {
			return nil, err
		}
	}
	return mono, nil
}

// sboxSumLazy computes each output bit as the ANF-weighted sum, splitting the terms into the
// factor monomials (degree-1, scale ~Delta) and the leaf monomials (degree-2, scale ~Delta^2).
// The two are summed separately: the leaf sum is relinearized+rescaled ONCE (deg-2 Delta^2 ->
// deg-1 Delta), then the two sums are combined. This defers all leaf relin/rescale to a single
// pair per output bit.
func (a *Evaluator) sboxSumLazy(mono map[int]*rlwe.Ciphertext, needed []int, coeffs [8][256]int, factors map[int]bool) (out ByteHE, err error) {
	monoCts := make([]*rlwe.Ciphertext, 0, len(needed))
	for _, S := range needed {
		monoCts = append(monoCts, mono[S])
	}
	utils.FlattenLevels(a.eval, monoCts)

	// addTerm accumulates c * mono[S] into *accPtr.
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
		var accF, accL *rlwe.Ciphertext // factors/bits (deg-1, ~Delta) and leaves (deg-2, ~Delta^2)
		for _, S := range needed {
			c := coeffs[i][S]
			if c == 0 {
				continue
			}
			if bits.OnesCount(uint(S)) < 2 || factors[S] {
				if err = addTerm(&accF, S, c); err != nil {
					return out, fmt.Errorf("sboxSumLazy bit %d factor S=%d: %w", i, S, err)
				}
			} else {
				if err = addTerm(&accL, S, c); err != nil {
					return out, fmt.Errorf("sboxSumLazy bit %d leaf S=%d: %w", i, S, err)
				}
			}
		}
		// Deferred relin + rescale of the leaf sum: one pair per output bit.
		if accL != nil {
			if err = a.eval.Relinearize(accL, accL); err != nil {
				return out, fmt.Errorf("sboxSumLazy bit %d relin: %w", i, err)
			}
			utils.Ops.Relin++
			if err = a.eval.Rescale(accL, accL); err != nil {
				return out, fmt.Errorf("sboxSumLazy bit %d rescale: %w", i, err)
			}
			utils.Ops.Rescale++
		}
		// Combine factor sum (deg-1, Delta) and leaf sum (deg-1, Delta after rescale).
		var acc *rlwe.Ciphertext
		switch {
		case accF != nil && accL != nil:
			utils.AlignLevels(a.eval, accF, accL)
			if err = a.eval.Add(accF, accL, accF); err != nil {
				return out, fmt.Errorf("sboxSumLazy bit %d combine: %w", i, err)
			}
			acc = accF
		case accF != nil:
			acc = accF
		case accL != nil:
			acc = accL
		default:
			return out, fmt.Errorf("sboxSumLazy bit %d: empty sum", i)
		}
		if a0 := coeffs[i][0]; a0 != 0 { // constant term a_empty
			if err = a.eval.Add(acc, a0, acc); err != nil {
				return out, fmt.Errorf("sboxSumLazy bit %d add const: %w", i, err)
			}
		}
		out[i] = acc
	}
	return out, nil
}

// SubByteV3 applies the AES S-box LAZILY: only the 61 monomials reused as a
// build factor are relinearized+rescaled (one pair each): the 186 leaf monomials stay degree-2
// at scale ~Delta^2 and their relin+rescale is deferred to one pair per output-bit accumulator.
// Cost: 61 + 8 = 69 relin and 69 rescale per byte, versus 247 for SubByteV2.
func (a *Evaluator) SubByteV3(inp ByteHE) (ByteHE, error) {
	coeffs := sboxCoeffs()
	needed := neededMonomials(coeffs)
	factors := factorMonomials(needed)
	mono, err := a.buildMonomialsLazy(inp, needed, factors)
	if err != nil {
		return ByteHE{}, err
	}
	return a.sboxSumLazy(mono, needed, coeffs, factors)
}

// buildMonomialsLazyV4 is like buildMonomialsLazy but keeps as lazy degree-2 products ONLY the
// leaf monomials of degree >= 4: the factors AND the low-degree (2,3) leaves are relinearized +
// rescaled at build time (utils.MulLeveled, degree-1).
func (a *Evaluator) buildMonomialsLazyV4(inp ByteHE, needed []int, factors map[int]bool) (map[int]*rlwe.Ciphertext, error) {
	mono := make(map[int]*rlwe.Ciphertext, len(needed)+8)
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
		if !factors[S] && bits.OnesCount(uint(S)) >= 4 {
			m, err = utils.MulLeveledLazy(a.eval, l, h) // deg>=4 leaf: defer relin + rescale (degree 2)
		} else {
			m, err = utils.MulLeveled(a.eval, l, h) // factor or deg<4 leaf: relin + rescale (degree 1)
		}
		if err != nil {
			return nil, fmt.Errorf("buildMonomialsLazyV4 S=%d: %w", S, err)
		}
		mono[S] = m
		return m, nil
	}
	for _, S := range needed {
		if _, err := build(S); err != nil {
			return nil, err
		}
	}
	return mono, nil
}

// sboxSumLazyV4 is like sboxSumLazy but the leaf accumulator collects ONLY the degree-2 leaves of
// degree >= 4 (the ones kept lazy by buildMonomialsLazyV4).
func (a *Evaluator) sboxSumLazyV4(mono map[int]*rlwe.Ciphertext, needed []int, coeffs [8][256]int, factors map[int]bool) (out ByteHE, err error) {
	monoCts := make([]*rlwe.Ciphertext, 0, len(needed))
	for _, S := range needed {
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
		var accF, accL *rlwe.Ciphertext // deg<4 or factor (deg-1, ~Delta) and deg>=4 leaves (deg-2, ~Delta^2)
		for _, S := range needed {
			c := coeffs[i][S]
			if c == 0 {
				continue
			}
			if !factors[S] && bits.OnesCount(uint(S)) >= 4 {
				if err = addTerm(&accL, S, c); err != nil {
					return out, fmt.Errorf("sboxSumLazyV4 bit %d leaf S=%d: %w", i, S, err)
				}
			} else {
				if err = addTerm(&accF, S, c); err != nil {
					return out, fmt.Errorf("sboxSumLazyV4 bit %d factor S=%d: %w", i, S, err)
				}
			}
		}
		if accL != nil {
			if err = a.eval.Relinearize(accL, accL); err != nil {
				return out, fmt.Errorf("sboxSumLazyV4 bit %d relin: %w", i, err)
			}
			utils.Ops.Relin++
			if err = a.eval.Rescale(accL, accL); err != nil {
				return out, fmt.Errorf("sboxSumLazyV4 bit %d rescale: %w", i, err)
			}
			utils.Ops.Rescale++
		}
		var acc *rlwe.Ciphertext
		switch {
		case accF != nil && accL != nil:
			utils.AlignLevels(a.eval, accF, accL)
			if err = a.eval.Add(accF, accL, accF); err != nil {
				return out, fmt.Errorf("sboxSumLazyV4 bit %d combine: %w", i, err)
			}
			acc = accF
		case accF != nil:
			acc = accF
		case accL != nil:
			acc = accL
		default:
			return out, fmt.Errorf("sboxSumLazyV4 bit %d: empty sum", i)
		}
		if a0 := coeffs[i][0]; a0 != 0 {
			if err = a.eval.Add(acc, a0, acc); err != nil {
				return out, fmt.Errorf("sboxSumLazyV4 bit %d add const: %w", i, err)
			}
		}
		out[i] = acc
	}
	return out, nil
}

// SubByteV4 applies the AES S-box with a less aggressive lazy strategy than SubByteV3: only the
// leaf monomials of degree >= 4 are kept lazy (degree-2, deferred): the factors and the
// degree-2/3 leaves are relinearized+rescaled at build time. Cost: 90 build + 8 output-bit
// accumulators = 98 relin and 98 rescale per byte.
func (a *Evaluator) SubByteV4(inp ByteHE) (ByteHE, error) {
	coeffs := sboxCoeffs()
	needed := neededMonomials(coeffs)
	factors := factorMonomials(needed)
	mono, err := a.buildMonomialsLazyV4(inp, needed, factors)
	if err != nil {
		return ByteHE{}, err
	}
	return a.sboxSumLazyV4(mono, needed, coeffs, factors)
}
