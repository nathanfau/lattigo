// Package utils holds small generic HE helpers shared across the nathanfau packages.
package utils

import (
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

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
