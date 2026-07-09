// Package utils holds small generic HE helpers shared across the nathanfau packages.
package utils

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// ============================================================================================
// QUICK REFERENCE
//
//   value? = what happens to the DECODED value
//     Y        preserved EXACTLY - DropLevel only (RNS limb truncation), adds no error
//     ~Y       preserved ALMOST exactly - a Rescale is involved, so its rounding noise is added
//     ~Y(bit)  ~Y but ONLY for p in {0,1} (p^2 = p); amplifies error near s=1, wrong otherwise
//     product  the func's job is to output x*y
//     combine  the func's job is to output a + i*b
//
//   scale? = what happens to the SCALE metadata
//     chain    stays on the canonical 2^Delta chain -> Add/Sub/Mul with other chain
//              ciphertexts is exact with NO manual scale matching
//     chain=   on the chain AND all outputs land on the SAME scale (ready to Add/Sub)
//     kept     original scale kept, may sit OFF the chain, so YOU must match it
//              (ForceScale) before combining with a chain ciphertext
//     Delta^2  degree-2 product scale (~Delta^2), not yet rescaled
//     set      forced to a scale target
//
//   level? = levels consumed by the call.
//
// 	 ops = relin + rescale actually performed.
//
//   func                  | value?  | scale?  | level?              					| ops (relin+resc)
//   ----------------------+---------+---------+----------------------------------------+-----------------
//   ForceScale(ct,s)      | ~Y      | set     | -1  (0 if already scaled at given s)	| 0+1
//   AlignLevels(a,b)      | Y       | kept    | min(a,b)       						| 0+0
//   FlattenLevels(cts)    | Y       | kept    | min(all)         						| 0+0
//   CombineReIm(a,b)      | combine | kept(a) | 0 (a and b are already align)			| 0+0
//   ----------------------+---------+---------+----------------------------------------+-----------------
//   SquareBit(ct)         | ~Y(bit) | chain   | -1                  					| 1+1
//   DescendBit(ct,l)      | ~Y(bit) | chain   | current_level - l = k 					| k*(1+1)
//   AlignBitLevels(a,b)   | ~Y(bit) | chain=  | min(a,b) = k   	    				| k*(1+1)
//   FlattenBitLevels(cts) | ~Y(bit) | chain=  | min(all) = k  	    					| k*(1+1)
//   ----------------------+---------+---------+----------------------------------------+-------------------
//   MulBits(x,y)          | product | chain   | -1 (+align squares) 					| 1+1 (+k*(1+1))
//   MulLeveled(a,b)       | product | kept    | -1                  					| 1+1
//   MulLeveledLazy(a,b)   | product | Delta^2 | 0  (output is deg-2)  					| 0+0
//
// Rule of thumb: a bit-sliced boolean circuit (homomorphic AES) uses the "chain" column
// (SquareBit / ... / MulBits) so every gate keeps scale = f(level) and Add/Sub need no
// fixups.
// Use the "kept" column (AlignLevels / MulLeveled / ...) when YOU manage scales
// yourself (RescaleTo / ForceScale to a fixed target).
// ============================================================================================

// OpCounter counts the expensive CKKS operations, used when designing differents SubByts versions.
type OpCounter struct {
	Relin   int
	Rescale int
}

var Ops OpCounter

func ResetOps() { Ops = OpCounter{} }

// ForceScale sets ct's scale to s, doing nothing if it is already exactly s.
func ForceScale(eval *ckks.Evaluator, ct *rlwe.Ciphertext, s rlwe.Scale) error {
	if ct.Scale.Cmp(s) == 0 {
		return nil
	}
	return eval.SetScale(ct, s)
}

// AlignLevels drops the level of whichever ciphertext is higher.
func AlignLevels(eval *ckks.Evaluator, a, b *rlwe.Ciphertext) {
	if a.Level() > b.Level() {
		eval.DropLevel(a, a.Level()-b.Level())
	} else if b.Level() > a.Level() {
		eval.DropLevel(b, b.Level()-a.Level())
	}
}

// CombineReIm returns a + i*b (levels aligned by DropLevel).
func CombineReIm(eval *ckks.Evaluator, a, b *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	ib := b.CopyNew()
	if err := eval.Mul(ib, complex(0, 1), ib); err != nil {
		return nil, fmt.Errorf("CombineReIm mul i: %w", err)
	}
	out := a.CopyNew()
	AlignLevels(eval, out, ib)
	if err := eval.Add(out, ib, out); err != nil {
		return nil, fmt.Errorf("CombineReIm add: %w", err)
	}
	return out, nil
}

// SquareBit computes ct <- ct^2 (MulRelin + Rescale), descending one level. For a bit p in
// {0,1}, p^2 = p, so the value is preserved while the scale follows the prime chain.
func SquareBit(eval *ckks.Evaluator, ct *rlwe.Ciphertext) error {
	if err := eval.MulRelin(ct, ct, ct); err != nil {
		return fmt.Errorf("SquareBit MulRelin: %w", err)
	}
	Ops.Relin++
	if err := eval.Rescale(ct, ct); err != nil {
		return fmt.Errorf("SquareBit Rescale: %w", err)
	}
	Ops.Rescale++
	return nil
}

// DescendBit lowers a bit ciphertext to level 'lvl' by repeated squaring (see SquareBit).
func DescendBit(eval *ckks.Evaluator, ct *rlwe.Ciphertext, lvl int) error {
	for ct.Level() > lvl {
		if err := SquareBit(eval, ct); err != nil {
			return fmt.Errorf("DescendBit: %w", err)
		}
	}
	return nil
}

// AlignBitLevels lowers whichever of the two bit ciphertexts is higher to the other's level
// by squaring. Unlike AlignLevels (DropLevel), the descended ciphertext keeps the canonical
// scale of its new level, so a and b end up at the same level AND the same scale.
func AlignBitLevels(eval *ckks.Evaluator, a, b *rlwe.Ciphertext) error {
	if a.Level() > b.Level() {
		return DescendBit(eval, a, b.Level())
	}
	if b.Level() > a.Level() {
		return DescendBit(eval, b, a.Level())
	}
	return nil
}

// FlattenBitLevels brings all bit ciphertexts to their common minimum level by squaring.
func FlattenBitLevels(eval *ckks.Evaluator, cts []*rlwe.Ciphertext) error {
	min := 1 << 30
	for _, ct := range cts {
		if l := ct.Level(); l < min {
			min = l
		}
	}
	for _, ct := range cts {
		if err := DescendBit(eval, ct, min); err != nil {
			return err
		}
	}
	return nil
}

// MulBits multiplies two bit ciphertexts, first aligning their levels by squaring (never
// DropLevel) so the product stays on the canonical scale chain, then MulRelin + Rescale.
//
// Lattigo's Mul already handles operands at different levels, but it descends the higher one
// by truncation (keeping its scale), which takes the product off the canonical scale chain
// and breaks later Add/Sub; squaring the higher operand (value-safe for a bit) avoids that.
func MulBits(eval *ckks.Evaluator, x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	p := x.CopyNew()
	q := y.CopyNew()
	if err := AlignBitLevels(eval, p, q); err != nil {
		return nil, fmt.Errorf("MulBits align: %w", err)
	}
	if err := eval.MulRelin(p, q, p); err != nil {
		return nil, fmt.Errorf("MulBits MulRelin: %w", err)
	}
	Ops.Relin++
	if err := eval.Rescale(p, p); err != nil {
		return nil, fmt.Errorf("MulBits Rescale: %w", err)
	}
	Ops.Rescale++
	return p, nil
}

// MulLeveled multiplies two ciphertexts the standard LEVELED way, WITHOUT squaring: the higher
// operand is brought to the lower operand's level by AlignLevels (DropLevel= no relin, no
// rescale, value-preserving), then a single MulRelin + Rescale.
func MulLeveled(eval *ckks.Evaluator, a, b *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	p := a.CopyNew()
	q := b.CopyNew()
	AlignLevels(eval, p, q)
	if err := eval.MulRelin(p, q, p); err != nil {
		return nil, fmt.Errorf("MulLeveled MulRelin: %w", err)
	}
	Ops.Relin++
	if err := eval.Rescale(p, p); err != nil {
		return nil, fmt.Errorf("MulLeveled Rescale: %w", err)
	}
	Ops.Rescale++
	return p, nil
}

// MulLeveledLazy multiplies two ciphertexts WITHOUT relinearizing or rescaling: it aligns their
// levels by DropLevel (AlignLevels) then Mul, returning a DEGREE-2 ciphertext at scale ~Delta^2.
func MulLeveledLazy(eval *ckks.Evaluator, a, b *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	p := a.CopyNew()
	q := b.CopyNew()
	AlignLevels(eval, p, q)
	out, err := eval.MulNew(p, q)
	if err != nil {
		return nil, fmt.Errorf("MulLeveledLazy Mul: %w", err)
	}
	return out, nil
}

// FlattenLevels brings all ciphertexts to their common minimum level using DropLevel only (no
// relin, no rescale, value-preserving).
func FlattenLevels(eval *ckks.Evaluator, cts []*rlwe.Ciphertext) {
	min := 1 << 30
	for _, ct := range cts {
		if l := ct.Level(); l < min {
			min = l
		}
	}
	for _, ct := range cts {
		if ct.Level() > min {
			eval.DropLevel(ct, ct.Level()-min)
		}
	}
}
