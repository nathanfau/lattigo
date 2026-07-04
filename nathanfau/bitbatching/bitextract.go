package bitbatching

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/circuits/common/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// BitExtract homomorphically extracts k bits from ct encrypting {omega^{m_s}}.
func BitExtract(params ckks.Parameters, eval *ckks.Evaluator, ct *rlwe.Ciphertext, k int) ([]*rlwe.Ciphertext, error) {
	if k < 1 {
		return nil, fmt.Errorf("BitExtract: k must be >= 1")
	}

	W := ct.Scale

	// keep the powers of 2 of X (ct) in memory so we don't recompute them every time
	pb := polynomial.NewPowerBasis(ct, bignum.Monomial)
	for i := 1; i < k; i++ {
		if err := pb.GenPower(1<<i, false, eval); err != nil {
			return nil, fmt.Errorf("BitExtract GenPower 2^%d: %w", i, err)
		}
	}

	// compute all the other powers
	for iter := k - 1; iter >= 1; iter-- {
		stride := 1 << (k - iter)
		L := 1 << (iter - 1)
		if L < 2 {
			continue
		}
		bsIdx, gsIdx := psPowerIndices(stride, L)
		for _, idx := range append(append([]int{}, bsIdx...), gsIdx...) {
			if idx <= 1 {
				continue
			}
			if err := pb.GenPower(idx, false, eval); err != nil {
				return nil, fmt.Errorf("BitExtract GenPower X^%d: %w", idx, err)
			}
		}
	}

	results := make([]*rlwe.Ciphertext, k)

	// exactly the BKSS algorithm
	for iter := k - 1; iter >= 1; iter-- {

		// computation of the polynomials P and Q (see polys.go)
		pkl := ComputePkl(k, iter)
		q := ComputeQkl(k, iter, pkl)

		prefix := pb.Value[1<<(k-iter-1)]

		var ctQ *rlwe.Ciphertext
		var err error

		if len(q) == 1 {
			ctQ = nil
		} else {
			stride := 1 << (k - iter)
			bsIdx, gsIdx := psPowerIndices(stride, len(q))
			if ctQ, err = evaluateCustomPS(eval, []complex128(q), bsIdx, gsIdx, &pb, W); err != nil {
				return nil, fmt.Errorf("BitExtract iter=%d evaluateCustomPS: %w", iter, err)
			}
		}

		var tmp *rlwe.Ciphertext
		if ctQ == nil {
			tmp = prefix.CopyNew()
			if err = eval.Mul(tmp, q[0], tmp); err != nil {
				return nil, fmt.Errorf("BitExtract iter=1 Mul: %w", err)
			}
			if err = eval.RescaleTo(tmp, W, tmp); err != nil {
				return nil, fmt.Errorf("BitExtract iter=1 RescaleTo: %w", err)
			}
		} else {
			p := prefix.CopyNew()
			qClone := ctQ.CopyNew()
			utils.AlignLevels(eval, p, qClone)
			tmp = p.CopyNew()
			if err = eval.MulRelin(p, qClone, tmp); err != nil {
				return nil, fmt.Errorf("BitExtract iter=%d MulRelin prefix×Q: %w", iter, err)
			}
			if err = eval.RescaleTo(tmp, W, tmp); err != nil {
				return nil, fmt.Errorf("BitExtract iter=%d RescaleTo prefix×Q: %w", iter, err)
			}
		}

		// recover the real part (the bit)
		if err = eval.Add(tmp, complex(0.25, 0), tmp); err != nil {
			return nil, fmt.Errorf("BitExtract iter=%d Add 1/4: %w", iter, err)
		}
		ctConj := tmp.CopyNew()
		if err = eval.Conjugate(ctConj, ctConj); err != nil {
			return nil, fmt.Errorf("BitExtract iter=%d Conjugate: %w", iter, err)
		}
		ctEll := tmp.CopyNew()
		if err = eval.Add(ctEll, ctConj, ctEll); err != nil {
			return nil, fmt.Errorf("BitExtract iter=%d Add+conj: %w", iter, err)
		}
		results[iter] = ctEll
	}

	ct0 := pb.Value[1<<(k-1)].CopyNew()
	if err := eval.Mul(ct0, complex(-0.5, 0), ct0); err != nil {
		return nil, fmt.Errorf("BitExtract LSB Mul: %w", err)
	}
	if err := eval.RescaleTo(ct0, W, ct0); err != nil {
		return nil, fmt.Errorf("BitExtract LSB RescaleTo: %w", err)
	}
	if err := eval.Add(ct0, complex(0.5, 0), ct0); err != nil {
		return nil, fmt.Errorf("BitExtract LSB Add: %w", err)
	}
	results[0] = ct0

	// normalize scale/level of every ct
	for i := range results {
		if err := utils.ForceScale(eval, results[i], W); err != nil {
			return nil, fmt.Errorf("BitExtract scale normalization bit %d: %w", i, err)
		}
	}
	minLvl := results[0].Level()
	for _, r := range results {
		if r.Level() < minLvl {
			minLvl = r.Level()
		}
	}
	for i := range results {
		if d := results[i].Level() - minLvl; d > 0 {
			eval.DropLevel(results[i], d)
		}
	}

	return results, nil
}

// babyStepM returns the number of baby-steps m (a power of 2) of the PS for a degree-L polynomial.
func babyStepM(L int) int {
	deg := L - 1
	if deg < 1 {
		return 1
	}
	logDegree := bits.Len(uint(deg))
	logSplit := bignum.OptimalSplit(logDegree)
	return 1 << logSplit
}

// psPowerIndices returns the baby-step and giant-step indices.
func psPowerIndices(stride, L int) (bs, gs []int) {
	m := babyStepM(L)
	maxBaby := m - 1
	if maxBaby > L-1 {
		maxBaby = L - 1
	}
	for i := 1; i <= maxBaby; i++ {
		bs = append(bs, stride*i)
	}
	for g := m; g <= L-1; g *= 2 {
		gs = append(gs, stride*g)
	}
	return
}

type psNode struct {
	ct  *rlwe.Ciphertext
	cst complex128
}

// evaluateCustomPS is a custom Paterson-Stockmeyer evaluation that reuses pb (the power
// basis)
// as noted above (l22), the power basis is not recomputed, which saves work.
/*


FLAG
this claim looks true but I should verify it



*/
func evaluateCustomPS(eval *ckks.Evaluator, coeffs []complex128, bsIdx, gsIdx []int, pb *polynomial.PowerBasis, W rlwe.Scale) (*rlwe.Ciphertext, error) {
	if len(bsIdx) == 0 {
		return nil, fmt.Errorf("evaluateCustomPS: empty bsIdx (insufficient degree)")
	}
	stride := bsIdx[0]
	var m int
	if len(gsIdx) > 0 {
		m = gsIdx[0] / stride
	} else {
		m = len(bsIdx) + 1
	}

	var rec func(c []complex128) (psNode, error)
	rec = func(c []complex128) (psNode, error) {
		n := len(c)

		if n <= m {
			return babyStepEval(eval, c, stride, pb, W)
		}

		g := m
		for g*2 <= n-1 {
			g *= 2
		}

		low, err := rec(c[:g])
		if err != nil {
			return psNode{}, err
		}
		high, err := rec(c[g:])
		if err != nil {
			return psNode{}, err
		}

		Zg := pb.Value[stride*g]
		if Zg == nil {
			return psNode{}, fmt.Errorf("evaluateCustomPS: giant X^%d missing from pb", stride*g)
		}
		var highCt *rlwe.Ciphertext
		if high.ct == nil {
			highCt = Zg.CopyNew()
			if err = eval.Mul(highCt, high.cst, highCt); err != nil {
				return psNode{}, fmt.Errorf("giant Mul const: %w", err)
			}
			if err = eval.RescaleTo(highCt, W, highCt); err != nil {
				return psNode{}, fmt.Errorf("giant RescaleTo const: %w", err)
			}
		} else {
			highCt = high.ct
			zc := Zg.CopyNew()
			utils.AlignLevels(eval, highCt, zc)
			if err = eval.MulRelin(highCt, zc, highCt); err != nil {
				return psNode{}, fmt.Errorf("giant MulRelin: %w", err)
			}
			if err = eval.RescaleTo(highCt, W, highCt); err != nil {
				return psNode{}, fmt.Errorf("giant RescaleTo: %w", err)
			}
		}

		if low.ct == nil {
			if err = eval.Add(highCt, low.cst, highCt); err != nil {
				return psNode{}, fmt.Errorf("giant Add low const: %w", err)
			}
			return psNode{ct: highCt}, nil
		}
		lc := low.ct
		utils.AlignLevels(eval, lc, highCt)
		if err = eval.Add(lc, highCt, lc); err != nil {
			return psNode{}, fmt.Errorf("giant Add low: %w", err)
		}
		return psNode{ct: lc}, nil
	}

	node, err := rec(coeffs)
	if err != nil {
		return nil, err
	}
	if node.ct == nil {
		return nil, fmt.Errorf("evaluateCustomPS: unexpected constant result")
	}
	return node.ct, nil
}

func babyStepEval(eval *ckks.Evaluator, c []complex128, stride int, pb *polynomial.PowerBasis, W rlwe.Scale) (psNode, error) {
	var acc *rlwe.Ciphertext
	for i := 1; i < len(c); i++ {
		if c[i] == 0 {
			continue
		}
		Yi := pb.Value[stride*i]
		if Yi == nil {
			return psNode{}, fmt.Errorf("babyStepEval: baby X^%d missing from pb", stride*i)
		}
		term := Yi.CopyNew()
		if err := eval.Mul(term, c[i], term); err != nil {
			return psNode{}, fmt.Errorf("baby Mul: %w", err)
		}
		if err := eval.RescaleTo(term, W, term); err != nil {
			return psNode{}, fmt.Errorf("baby RescaleTo: %w", err)
		}
		if acc == nil {
			acc = term
		} else {
			utils.AlignLevels(eval, acc, term)
			if err := eval.Add(acc, term, acc); err != nil {
				return psNode{}, fmt.Errorf("baby Add: %w", err)
			}
		}
	}
	if acc == nil {
		return psNode{ct: nil, cst: c[0]}, nil
	}
	if c[0] != 0 {
		if err := eval.Add(acc, c[0], acc); err != nil {
			return psNode{}, fmt.Errorf("baby Add const: %w", err)
		}
	}
	return psNode{ct: acc}, nil
}
