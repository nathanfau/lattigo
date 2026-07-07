// Package utils holds small generic HE helpers shared across the nathanfau packages.
package utils

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

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

// The functions below descend a ciphertext to a lower level (or align two/many to a common
// level) by SQUARING, as an alternative to AlignLevels/DropLevel above. The two families
// solve different problems and are NOT interchangeable:
//
//   - AlignLevels (DropLevel) is VALUE-preserving, but keeps the ORIGINAL SCALE.
//     The dropped ciphertext no longer sits on the "scale reached by successive rescales"
//     chain, so it can only be combined (Add/Sub) with a ciphertext whose scale has been
//     matched by hand (see ForceScale), otherwise the scales mismatch.
//
//   - SquareBit, ...: descend by squarring (MulRelin + Rescale). This is VALUE-preserving only for
//     a bit p in {0,1}, but the Rescale keeps the scale on the canonical
//     chain: the scale at a given level is exactly what repeated squaring produces. Two bit
//     ciphertexts brought to the same level this way share the SAME SCALE, so Add/Sub/Mul
//     between them are exact with no ForceScale. This is what lets a bit-sliced boolean
//     circuit (e.g. the homomorphic AES) keep the invariant scale = f(level) across gates.
//
// Rule : Use AlignLevels when you MANAGE SCALES YOURSELF (RescaleTo/ForceScale to a
// fixed target) and use the squaring variants when your ciphertexts encrypt bits and you want each level
// drop to also fix up the scale automatically (one operation does both).

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
