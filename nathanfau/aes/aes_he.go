package aes

import (
	"fmt"
	"sort"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// ByteHE is a bit-sliced encrypted byte: 8 ciphertexts, one per bit.
type ByteHE = [8]*rlwe.Ciphertext

// Evaluator runs the bit-sliced AES round functions over CKKS ciphertexts.
type Evaluator struct {
	eval *ckks.Evaluator
}

func NewEvaluator(eval *ckks.Evaluator) *Evaluator {
	return &Evaluator{eval: eval}
}

// Xor computes x XOR y = (x - y)^2 between bits (ciphertexts in {0,1})
func (a *Evaluator) Xor(x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	xx := x.CopyNew()
	yy := y.CopyNew()

	utils.AlignLevels(a.eval, xx, yy)
	if err := a.eval.Sub(xx, yy, xx); err != nil { // x - y
		return nil, fmt.Errorf("xor x-y: %w", err)
	}
	if err := a.eval.MulRelin(xx, xx, xx); err != nil { // (x - y)^2
		return nil, fmt.Errorf("xor square: %w", err)
	}
	if err := a.eval.Rescale(xx, xx); err != nil {
		return nil, fmt.Errorf("xor rescale: %w", err)
	}
	return xx, nil
}

type mcLeaf struct{ off, bit int }

// mixColLeaves generates, for each output D_d and each bit j, the reduced-parity (sorted for
// reproducibility) list of input leaves. In encryption form:
//
//	D_d = xtime(b_p ^ b_q) ^ b_e0 ^ b_e1 ^ b_e2
//
// only kept for reproducibility.
func mixColLeaves() [4][8][]mcLeaf {
	xmap := [8][]int{{7}, {0, 7}, {1}, {2, 7}, {3, 7}, {4}, {5}, {6}}
	type ddef struct {
		p, q   int
		extras [3]int
	}
	dd := [4]ddef{
		{0, 1, [3]int{1, 2, 3}},
		{1, 2, [3]int{2, 3, 0}},
		{2, 3, [3]int{3, 0, 1}},
		{3, 0, [3]int{0, 1, 2}},
	}
	var res [4][8][]mcLeaf
	for d := 0; d < 4; d++ {
		for j := 0; j < 8; j++ {
			cnt := map[mcLeaf]int{}
			for _, k := range xmap[j] {
				cnt[mcLeaf{dd[d].p, k}]++
				cnt[mcLeaf{dd[d].q, k}]++
			}
			for _, m := range dd[d].extras {
				cnt[mcLeaf{m, j}]++
			}
			for lf, c := range cnt {
				if c%2 == 1 {
					res[d][j] = append(res[d][j], lf)
				}
			}
			sort.Slice(res[d][j], func(x, y int) bool {
				a, b := res[d][j][x], res[d][j][y]
				if a.off != b.off {
					return a.off < b.off
				}
				return a.bit < b.bit
			})
		}
	}
	return res
}

// mcLeavesTable is the hard-coded output of mixColLeaves
var mcLeavesTable = [4][8][]mcLeaf{
	{{{0, 7}, {1, 0}, {1, 7}, {2, 0}, {3, 0}}, {{0, 0}, {0, 7}, {1, 0}, {1, 1}, {1, 7}, {2, 1}, {3, 1}}, {{0, 1}, {1, 1}, {1, 2}, {2, 2}, {3, 2}}, {{0, 2}, {0, 7}, {1, 2}, {1, 3}, {1, 7}, {2, 3}, {3, 3}}, {{0, 3}, {0, 7}, {1, 3}, {1, 4}, {1, 7}, {2, 4}, {3, 4}}, {{0, 4}, {1, 4}, {1, 5}, {2, 5}, {3, 5}}, {{0, 5}, {1, 5}, {1, 6}, {2, 6}, {3, 6}}, {{0, 6}, {1, 6}, {1, 7}, {2, 7}, {3, 7}}}, // D0
	{{{0, 0}, {1, 7}, {2, 0}, {2, 7}, {3, 0}}, {{0, 1}, {1, 0}, {1, 7}, {2, 0}, {2, 1}, {2, 7}, {3, 1}}, {{0, 2}, {1, 1}, {2, 1}, {2, 2}, {3, 2}}, {{0, 3}, {1, 2}, {1, 7}, {2, 2}, {2, 3}, {2, 7}, {3, 3}}, {{0, 4}, {1, 3}, {1, 7}, {2, 3}, {2, 4}, {2, 7}, {3, 4}}, {{0, 5}, {1, 4}, {2, 4}, {2, 5}, {3, 5}}, {{0, 6}, {1, 5}, {2, 5}, {2, 6}, {3, 6}}, {{0, 7}, {1, 6}, {2, 6}, {2, 7}, {3, 7}}}, // D1
	{{{0, 0}, {1, 0}, {2, 7}, {3, 0}, {3, 7}}, {{0, 1}, {1, 1}, {2, 0}, {2, 7}, {3, 0}, {3, 1}, {3, 7}}, {{0, 2}, {1, 2}, {2, 1}, {3, 1}, {3, 2}}, {{0, 3}, {1, 3}, {2, 2}, {2, 7}, {3, 2}, {3, 3}, {3, 7}}, {{0, 4}, {1, 4}, {2, 3}, {2, 7}, {3, 3}, {3, 4}, {3, 7}}, {{0, 5}, {1, 5}, {2, 4}, {3, 4}, {3, 5}}, {{0, 6}, {1, 6}, {2, 5}, {3, 5}, {3, 6}}, {{0, 7}, {1, 7}, {2, 6}, {3, 6}, {3, 7}}}, // D2
	{{{0, 0}, {0, 7}, {1, 0}, {2, 0}, {3, 7}}, {{0, 0}, {0, 1}, {0, 7}, {1, 1}, {2, 1}, {3, 0}, {3, 7}}, {{0, 1}, {0, 2}, {1, 2}, {2, 2}, {3, 1}}, {{0, 2}, {0, 3}, {0, 7}, {1, 3}, {2, 3}, {3, 2}, {3, 7}}, {{0, 3}, {0, 4}, {0, 7}, {1, 4}, {2, 4}, {3, 3}, {3, 7}}, {{0, 4}, {0, 5}, {1, 5}, {2, 5}, {3, 4}}, {{0, 5}, {0, 6}, {1, 6}, {2, 6}, {3, 5}}, {{0, 6}, {0, 7}, {1, 7}, {2, 7}, {3, 6}}}, // D3
}

// xorTree XORs a list of ciphertexts via a balanced tree.
func (a *Evaluator) xorTree(cts []*rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return utils.ReduceBalanced(cts, a.Xor)
}

var (
	shiftRowsSrc  = [8]int{0, 5, 2, 7, 4, 1, 6, 3}
	shiftRowsSwap = [8]bool{false, false, true, true, false, true, true, false}
)

// ShiftRows applies AES ShiftRows at the Algo1 pause
func (a *Evaluator) ShiftRows(reals, imags [16]*rlwe.Ciphertext) (re, im [16]*rlwe.Ciphertext) {
	for j := 0; j < 8; j++ {
		g := shiftRowsSrc[j]
		if shiftRowsSwap[j] {
			re[2*j], re[2*j+1] = imags[2*g], imags[2*g+1]
			im[2*j], im[2*j+1] = reals[2*g], reals[2*g+1]
		} else {
			re[2*j], re[2*j+1] = reals[2*g], reals[2*g+1]
			im[2*j], im[2*j+1] = imags[2*g], imags[2*g+1]
		}
	}
	return re, im
}

// MixColumns applies AES MixColumns to a blockpack's way packed group of ct
func (a *Evaluator) MixColumns(p [8]ByteHE) ([8]ByteHE, error) {
	var out [8]ByteHE
	for quad := 0; quad < 2; quad++ {
		base := 4 * quad
		for d := 0; d < 4; d++ {
			var ob ByteHE
			for j := 0; j < 8; j++ {
				leaves := mcLeavesTable[d][j]
				cts := make([]*rlwe.Ciphertext, len(leaves))
				for n, lf := range leaves {
					cts[n] = p[base+lf.off][lf.bit]
				}
				r, err := a.xorTree(cts)
				if err != nil {
					return out, fmt.Errorf("MixColumns quad %d D%d bit %d: %w", quad, d, j, err)
				}
				ob[j] = r
			}
			out[base+d] = ob
		}
	}
	return out, nil
}
