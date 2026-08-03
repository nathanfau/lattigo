package cleaning

// The composed routes of xorclean.go, written out from A to Z: no Cleaning, no aes XOR, no utils
// helper -- only the primitive evaluator operations, so nothing hides behind a function call.
// h(t) = 3t^2 - 2t^3, evaluated as t^2 * (3 - 2t). Three levels each.
//
// Only Cleaning is covered: the two smoother polynomials would need their own expansion.

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func cleanThenXorSqByHand(_ ckks.Parameters, eval *ckks.Evaluator, x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	// h(x)
	hx := x.CopyNew()
	if err := eval.MulRelin(hx, hx, hx); err != nil { // x^2
		return nil, fmt.Errorf("byHand x^2: %w", err)
	}
	if err := eval.Rescale(hx, hx); err != nil {
		return nil, fmt.Errorf("byHand x^2 rescale: %w", err)
	}
	lx := x.CopyNew()
	if err := eval.Mul(lx, -2, lx); err != nil { // 3 - 2x
		return nil, fmt.Errorf("byHand -2x: %w", err)
	}
	if err := eval.Add(lx, 3, lx); err != nil {
		return nil, fmt.Errorf("byHand 3-2x: %w", err)
	}
	if lx.Level() > hx.Level() {
		eval.DropLevel(lx, lx.Level()-hx.Level())
	}
	if err := eval.MulRelin(hx, lx, hx); err != nil {
		return nil, fmt.Errorf("byHand h(x): %w", err)
	}
	if err := eval.Rescale(hx, hx); err != nil {
		return nil, fmt.Errorf("byHand h(x) rescale: %w", err)
	}

	// h(y)
	hy := y.CopyNew()
	if err := eval.MulRelin(hy, hy, hy); err != nil { // y^2
		return nil, fmt.Errorf("byHand y^2: %w", err)
	}
	if err := eval.Rescale(hy, hy); err != nil {
		return nil, fmt.Errorf("byHand y^2 rescale: %w", err)
	}
	ly := y.CopyNew()
	if err := eval.Mul(ly, -2, ly); err != nil { // 3 - 2y
		return nil, fmt.Errorf("byHand -2y: %w", err)
	}
	if err := eval.Add(ly, 3, ly); err != nil {
		return nil, fmt.Errorf("byHand 3-2y: %w", err)
	}
	if ly.Level() > hy.Level() {
		eval.DropLevel(ly, ly.Level()-hy.Level())
	}
	if err := eval.MulRelin(hy, ly, hy); err != nil {
		return nil, fmt.Errorf("byHand h(y): %w", err)
	}
	if err := eval.Rescale(hy, hy); err != nil {
		return nil, fmt.Errorf("byHand h(y) rescale: %w", err)
	}

	// (h(x) - h(y))^2
	if hx.Level() > hy.Level() {
		eval.DropLevel(hx, hx.Level()-hy.Level())
	} else if hy.Level() > hx.Level() {
		eval.DropLevel(hy, hy.Level()-hx.Level())
	}
	if err := eval.Sub(hx, hy, hx); err != nil {
		return nil, fmt.Errorf("byHand sub: %w", err)
	}
	if err := eval.MulRelin(hx, hx, hx); err != nil {
		return nil, fmt.Errorf("byHand square: %w", err)
	}
	return hx, eval.Rescale(hx, hx)
}

func cleanThenXorNoSqByHand(_ ckks.Parameters, eval *ckks.Evaluator, x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	// h(x)
	hx := x.CopyNew()
	if err := eval.MulRelin(hx, hx, hx); err != nil {
		return nil, fmt.Errorf("byHand x^2: %w", err)
	}
	if err := eval.Rescale(hx, hx); err != nil {
		return nil, fmt.Errorf("byHand x^2 rescale: %w", err)
	}
	lx := x.CopyNew()
	if err := eval.Mul(lx, -2, lx); err != nil {
		return nil, fmt.Errorf("byHand -2x: %w", err)
	}
	if err := eval.Add(lx, 3, lx); err != nil {
		return nil, fmt.Errorf("byHand 3-2x: %w", err)
	}
	if lx.Level() > hx.Level() {
		eval.DropLevel(lx, lx.Level()-hx.Level())
	}
	if err := eval.MulRelin(hx, lx, hx); err != nil {
		return nil, fmt.Errorf("byHand h(x): %w", err)
	}
	if err := eval.Rescale(hx, hx); err != nil {
		return nil, fmt.Errorf("byHand h(x) rescale: %w", err)
	}

	// h(y)
	hy := y.CopyNew()
	if err := eval.MulRelin(hy, hy, hy); err != nil {
		return nil, fmt.Errorf("byHand y^2: %w", err)
	}
	if err := eval.Rescale(hy, hy); err != nil {
		return nil, fmt.Errorf("byHand y^2 rescale: %w", err)
	}
	ly := y.CopyNew()
	if err := eval.Mul(ly, -2, ly); err != nil {
		return nil, fmt.Errorf("byHand -2y: %w", err)
	}
	if err := eval.Add(ly, 3, ly); err != nil {
		return nil, fmt.Errorf("byHand 3-2y: %w", err)
	}
	if ly.Level() > hy.Level() {
		eval.DropLevel(ly, ly.Level()-hy.Level())
	}
	if err := eval.MulRelin(hy, ly, hy); err != nil {
		return nil, fmt.Errorf("byHand h(y): %w", err)
	}
	if err := eval.Rescale(hy, hy); err != nil {
		return nil, fmt.Errorf("byHand h(y) rescale: %w", err)
	}

	// 1/2 - (1/2 - h(x))(1 - 2h(y))
	if hx.Level() > hy.Level() {
		eval.DropLevel(hx, hx.Level()-hy.Level())
	} else if hy.Level() > hx.Level() {
		eval.DropLevel(hy, hy.Level()-hx.Level())
	}
	if err := eval.Mul(hx, -1, hx); err != nil {
		return nil, fmt.Errorf("byHand -h(x): %w", err)
	}
	if err := eval.Add(hx, 0.5, hx); err != nil {
		return nil, fmt.Errorf("byHand 1/2-h(x): %w", err)
	}
	if err := eval.Mul(hy, -2, hy); err != nil {
		return nil, fmt.Errorf("byHand -2h(y): %w", err)
	}
	if err := eval.Add(hy, 1, hy); err != nil {
		return nil, fmt.Errorf("byHand 1-2h(y): %w", err)
	}
	if err := eval.MulRelin(hx, hy, hx); err != nil {
		return nil, fmt.Errorf("byHand mul: %w", err)
	}
	if err := eval.Rescale(hx, hx); err != nil {
		return nil, fmt.Errorf("byHand rescale: %w", err)
	}
	if err := eval.Mul(hx, -1, hx); err != nil {
		return nil, fmt.Errorf("byHand negate: %w", err)
	}
	if err := eval.Add(hx, 0.5, hx); err != nil {
		return nil, fmt.Errorf("byHand 1/2-uv: %w", err)
	}
	return hx, nil
}

func cleanOneThenXorSqByHand(_ ckks.Parameters, eval *ckks.Evaluator, x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	// h(x)
	hx := x.CopyNew()
	if err := eval.MulRelin(hx, hx, hx); err != nil {
		return nil, fmt.Errorf("byHand x^2: %w", err)
	}
	if err := eval.Rescale(hx, hx); err != nil {
		return nil, fmt.Errorf("byHand x^2 rescale: %w", err)
	}
	lx := x.CopyNew()
	if err := eval.Mul(lx, -2, lx); err != nil {
		return nil, fmt.Errorf("byHand -2x: %w", err)
	}
	if err := eval.Add(lx, 3, lx); err != nil {
		return nil, fmt.Errorf("byHand 3-2x: %w", err)
	}
	if lx.Level() > hx.Level() {
		eval.DropLevel(lx, lx.Level()-hx.Level())
	}
	if err := eval.MulRelin(hx, lx, hx); err != nil {
		return nil, fmt.Errorf("byHand h(x): %w", err)
	}
	if err := eval.Rescale(hx, hx); err != nil {
		return nil, fmt.Errorf("byHand h(x) rescale: %w", err)
	}

	// (h(x) - y)^2, y untouched
	yy := y.CopyNew()
	if hx.Level() > yy.Level() {
		eval.DropLevel(hx, hx.Level()-yy.Level())
	} else if yy.Level() > hx.Level() {
		eval.DropLevel(yy, yy.Level()-hx.Level())
	}
	if err := eval.Sub(hx, yy, hx); err != nil {
		return nil, fmt.Errorf("byHand sub: %w", err)
	}
	if err := eval.MulRelin(hx, hx, hx); err != nil {
		return nil, fmt.Errorf("byHand square: %w", err)
	}
	return hx, eval.Rescale(hx, hx)
}

func cleanOneThenXorNoSqByHand(_ ckks.Parameters, eval *ckks.Evaluator, x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	// h(x)
	hx := x.CopyNew()
	if err := eval.MulRelin(hx, hx, hx); err != nil {
		return nil, fmt.Errorf("byHand x^2: %w", err)
	}
	if err := eval.Rescale(hx, hx); err != nil {
		return nil, fmt.Errorf("byHand x^2 rescale: %w", err)
	}
	lx := x.CopyNew()
	if err := eval.Mul(lx, -2, lx); err != nil {
		return nil, fmt.Errorf("byHand -2x: %w", err)
	}
	if err := eval.Add(lx, 3, lx); err != nil {
		return nil, fmt.Errorf("byHand 3-2x: %w", err)
	}
	if lx.Level() > hx.Level() {
		eval.DropLevel(lx, lx.Level()-hx.Level())
	}
	if err := eval.MulRelin(hx, lx, hx); err != nil {
		return nil, fmt.Errorf("byHand h(x): %w", err)
	}
	if err := eval.Rescale(hx, hx); err != nil {
		return nil, fmt.Errorf("byHand h(x) rescale: %w", err)
	}

	// 1/2 - (1/2 - h(x))(1 - 2y), y untouched
	yy := y.CopyNew()
	if hx.Level() > yy.Level() {
		eval.DropLevel(hx, hx.Level()-yy.Level())
	} else if yy.Level() > hx.Level() {
		eval.DropLevel(yy, yy.Level()-hx.Level())
	}
	if err := eval.Mul(hx, -1, hx); err != nil {
		return nil, fmt.Errorf("byHand -h(x): %w", err)
	}
	if err := eval.Add(hx, 0.5, hx); err != nil {
		return nil, fmt.Errorf("byHand 1/2-h(x): %w", err)
	}
	if err := eval.Mul(yy, -2, yy); err != nil {
		return nil, fmt.Errorf("byHand -2y: %w", err)
	}
	if err := eval.Add(yy, 1, yy); err != nil {
		return nil, fmt.Errorf("byHand 1-2y: %w", err)
	}
	if err := eval.MulRelin(hx, yy, hx); err != nil {
		return nil, fmt.Errorf("byHand mul: %w", err)
	}
	if err := eval.Rescale(hx, hx); err != nil {
		return nil, fmt.Errorf("byHand rescale: %w", err)
	}
	if err := eval.Mul(hx, -1, hx); err != nil {
		return nil, fmt.Errorf("byHand negate: %w", err)
	}
	if err := eval.Add(hx, 0.5, hx); err != nil {
		return nil, fmt.Errorf("byHand 1/2-uv: %w", err)
	}
	return hx, nil
}
