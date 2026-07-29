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

// MulLeveled multiplies two ciphertexts the standard LEVELED way: the higher operand is brought
// to the lower operand's level by AlignLevels (DropLevel = no relin, no rescale, value-preserving),
// then a single MulRelin + Rescale.
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

// ReduceBalanced folds a list into a single element with a BALANCED tree, i.e. depth
// ceil(log2 n) instead of the n-1 of a left fold, which is what keeps the XOR trees shallow.
func ReduceBalanced[T any](items []T, combine func(T, T) (T, error)) (T, error) {
	var zero T
	if len(items) == 0 {
		return zero, fmt.Errorf("ReduceBalanced: empty list")
	}
	cur := make([]T, len(items))
	copy(cur, items)
	for len(cur) > 1 {
		var next []T
		for i := 0; i+1 < len(cur); i += 2 {
			x, err := combine(cur[i], cur[i+1])
			if err != nil {
				return zero, err
			}
			next = append(next, x)
		}
		if len(cur)%2 == 1 {
			next = append(next, cur[len(cur)-1])
		}
		cur = next
	}
	return cur[0], nil
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
