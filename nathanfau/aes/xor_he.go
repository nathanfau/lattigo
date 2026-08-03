package aes

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
)

// XorSq computes x XOR y = (x - y)^2
func (a *Evaluator) XorSq(x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	xx := x.CopyNew()
	yy := y.CopyNew()

	utils.AlignLevels(a.eval, xx, yy)
	if err := a.eval.Sub(xx, yy, xx); err != nil { // x - y
		return nil, fmt.Errorf("XorSq x-y: %w", err)
	}
	if err := a.eval.MulRelin(xx, xx, xx); err != nil { // (x - y)^2
		return nil, fmt.Errorf("XorSq square: %w", err)
	}
	if err := a.eval.Rescale(xx, xx); err != nil {
		return nil, fmt.Errorf("XorSq rescale: %w", err)
	}
	return xx, nil
}

// XorNoSq computes x XOR y = x + y - 2xy, evaluated as 1/2 - (1/2 - x)(1 - 2y)
func (a *Evaluator) XorNoSq(x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	u := x.CopyNew()
	v := y.CopyNew()
	utils.AlignLevels(a.eval, u, v)

	if err := a.eval.Mul(u, -1, u); err != nil { // u = 1/2 - x
		return nil, fmt.Errorf("XorNoSq -x: %w", err)
	}
	if err := a.eval.Add(u, 0.5, u); err != nil {
		return nil, fmt.Errorf("XorNoSq 1/2-x: %w", err)
	}
	if err := a.eval.Mul(v, -2, v); err != nil { // v = 1 - 2y
		return nil, fmt.Errorf("XorNoSq -2y: %w", err)
	}
	if err := a.eval.Add(v, 1, v); err != nil {
		return nil, fmt.Errorf("XorNoSq 1-2y: %w", err)
	}
	if err := a.eval.MulRelin(u, v, u); err != nil { // uv
		return nil, fmt.Errorf("XorNoSq mul: %w", err)
	}
	if err := a.eval.Rescale(u, u); err != nil {
		return nil, fmt.Errorf("XorNoSq rescale: %w", err)
	}
	if err := a.eval.Mul(u, -1, u); err != nil { // 1/2 - uv
		return nil, fmt.Errorf("XorNoSq negate: %w", err)
	}
	if err := a.eval.Add(u, 0.5, u); err != nil {
		return nil, fmt.Errorf("XorNoSq 1/2-uv: %w", err)
	}
	return u, nil
}

// XorSqPlain is XorSq for ct*pt
func (a *Evaluator) XorSqPlain(x *rlwe.Ciphertext, bits []float64) (*rlwe.Ciphertext, error) {
	out := x.CopyNew()
	if err := a.eval.Sub(out, bits, out); err != nil {
		return nil, fmt.Errorf("XorSqPlain sub: %w", err)
	}
	if err := a.eval.MulRelin(out, out, out); err != nil {
		return nil, fmt.Errorf("XorSqPlain square: %w", err)
	}
	if err := a.eval.Rescale(out, out); err != nil {
		return nil, fmt.Errorf("XorSqPlain rescale: %w", err)
	}
	return out, nil
}

// XorNoSqPlain is XorNoSq for ct*pt
func (a *Evaluator) XorNoSqPlain(x *rlwe.Ciphertext, bits []float64) (*rlwe.Ciphertext, error) {
	m := make([]float64, len(bits))
	for i, b := range bits {
		m[i] = 1 - 2*b
	}
	out := x.CopyNew()
	if err := a.eval.Mul(out, m, out); err != nil {
		return nil, fmt.Errorf("XorNoSqPlain mul: %w", err)
	}
	if err := a.eval.Rescale(out, out); err != nil {
		return nil, fmt.Errorf("XorNoSqPlain rescale: %w", err)
	}
	if err := a.eval.Add(out, bits, out); err != nil {
		return nil, fmt.Errorf("XorNoSqPlain add: %w", err)
	}
	return out, nil
}
