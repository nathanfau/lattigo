package utils

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// RescalingFactor is what Rescale divides a ciphertext at this level by: the product of the primes
// it consumes, which is one of them unless the parameters rescale over several.
func RescalingFactor(params ckks.Parameters, level int) rlwe.Scale {
	ringQ := params.RingQ()
	s := rlwe.NewScale(1)
	for i := 0; i < params.LevelsConsumedPerRescaling(); i++ {
		s = s.Mul(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
	}
	return s
}

// MulRescaleTo multiplies ct by the real constant c and rescales it, landing on target rather than
// on ct's own scale: a scale hop that costs nothing beyond the Rescale. The constant is folded into
// the integer M = round(|c|.target.q/ct.Scale) and the declared scale M/|c|, so the value is
// exactly c times the input, and the scale lands on target up to 1/M.
func MulRescaleTo(eval *ckks.Evaluator, ct *rlwe.Ciphertext, c float64, target rlwe.Scale) error {
	if c == 0 {
		return fmt.Errorf("MulRescaleTo: zero constant")
	}
	absC := new(big.Float).SetPrec(rlwe.ScalePrecision).SetFloat64(c)
	absC.Abs(absC)

	d := target.Mul(RescalingFactor(*eval.GetParameters(), ct.Level())).Div(ct.Scale)
	mf := new(big.Float).SetPrec(rlwe.ScalePrecision).Mul(absC, &d.Value)
	mf.Add(mf, big.NewFloat(0.5))
	m, _ := mf.Int(nil)
	if m.Sign() == 0 {
		return fmt.Errorf("MulRescaleTo: |%g| at scale 2^%.2f rounds to 0", c, d.Log2())
	}
	decl := rlwe.NewScale(new(big.Float).SetPrec(rlwe.ScalePrecision).Quo(new(big.Float).SetInt(m), absC))
	if c < 0 {
		m.Neg(m)
	}

	if err := eval.Mul(ct, m, ct); err != nil { // an integer: the scale does not move
		return fmt.Errorf("MulRescaleTo Mul: %w", err)
	}
	ct.Scale = ct.Scale.Mul(decl)
	if err := eval.Rescale(ct, ct); err != nil {
		return fmt.Errorf("MulRescaleTo Rescale: %w", err)
	}
	return nil
}
