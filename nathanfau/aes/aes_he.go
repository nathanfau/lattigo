package aes

import (
	"fmt"
	"sort"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// ByteHE is a bit-sliced encrypted byte: 8 ciphertexts, one per bit.
type ByteHE [8]*rlwe.Ciphertext

// StateHE is a bit-sliced encrypted AES state: 16 bytes.
type StateHE [16]ByteHE

// Evaluator runs the bit-sliced AES round functions over CKKS ciphertexts.
//
// When noTrick is set, the XOR-based ops (xor, and everything built on it: xorByte, AddRoundKey,
// xtime, MixColumns) drop the squaring/scale-chain "trick" and use the standard leveled path
// (DropLevel alignment + MulRelin/Rescale) instead.
type Evaluator struct {
	eval    *ckks.Evaluator
	noTrick bool
}

// NewEvaluator builds the ~optimized (trick) evaluator
func NewEvaluator(eval *ckks.Evaluator) *Evaluator {
	return &Evaluator{eval: eval}
}

// NewEvaluatorNoTrick builds the ~standard evaluator: the XOR-based ops use the standard
// leveled path
func NewEvaluatorNoTrick(eval *ckks.Evaluator) *Evaluator {
	return &Evaluator{eval: eval, noTrick: true}
}

// xor computes x XOR y between bits (ciphertexts in {0,1}). In trick mode it is (x - y)^2 with
// level alignment by squaring (utils.AlignBitLevels + utils.SquareBit)
func (a *Evaluator) xor(x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	if a.noTrick {
		return a.xorStd(x, y)
	}

	xx := x.CopyNew()
	yy := y.CopyNew()

	if err := utils.AlignBitLevels(a.eval, xx, yy); err != nil {
		return nil, fmt.Errorf("xor align: %w", err)
	}
	if err := a.eval.Sub(xx, yy, xx); err != nil {
		return nil, fmt.Errorf("xx-yy: %w", err)
	}
	if err := utils.SquareBit(a.eval, xx); err != nil {
		return nil, fmt.Errorf("xor square: %w", err)
	}
	return xx, nil
}

// xorStd computes x XOR y = (x - y)^2 for bits: the SAME square as the trick xor, but the standard
// leveled way. Levels are aligned by DropLevel (utils.AlignLevels, no squaring) instead of the
// trick's AlignBitLevels, and the difference is squared with a plain MulRelin + Rescale instead of
// utils.SquareBit. Consumes one level, like the trick variant, but the scale is not kept on the
// canonical chain.
func (a *Evaluator) xorStd(x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	xx := x.CopyNew()
	yy := y.CopyNew()

	utils.AlignLevels(a.eval, xx, yy)
	if err := a.eval.Sub(xx, yy, xx); err != nil { // x - y
		return nil, fmt.Errorf("xorStd x-y: %w", err)
	}
	if err := a.eval.MulRelin(xx, xx, xx); err != nil { // (x - y)^2
		return nil, fmt.Errorf("xorStd square: %w", err)
	}
	if err := a.eval.Rescale(xx, xx); err != nil {
		return nil, fmt.Errorf("xorStd rescale: %w", err)
	}
	return xx, nil
}

// xorByte XORs two bytes bit by bit.
func (a *Evaluator) xorByte(x, y ByteHE) (ByteHE, error) {
	var out ByteHE
	for i := 0; i < 8; i++ {
		c, err := a.xor(x[i], y[i])
		if err != nil {
			return out, fmt.Errorf("xorByte bit %d: %w", i, err)
		}
		out[i] = c
	}
	return out, nil
}

// reduceBalanced folds a list into a single element with a BALANCED tree = depth ceil(log2 n)
func reduceBalanced[T any](items []T, combine func(T, T) (T, error)) (T, error) {
	var zero T
	if len(items) == 0 {
		return zero, fmt.Errorf("reduceBalanced: empty list")
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

// xorBytes XORs k>2 bytes via a balanced tree.
func (a *Evaluator) xorBytes(bytes ...ByteHE) (ByteHE, error) {
	return reduceBalanced(bytes, a.xorByte)
}

// xtime multiplies by X in GF(256):
//
//	out0=a7  out1=a0^a7  out2=a1  out3=a2^a7  out4=a3^a7  out5=a4  out6=a5  out7=a6
func (a *Evaluator) xtime(x ByteHE) (ByteHE, error) {
	var out ByteHE
	var err error
	out[0] = x[7].CopyNew()
	if out[1], err = a.xor(x[0], x[7]); err != nil {
		return out, fmt.Errorf("xtime out1: %w", err)
	}
	out[2] = x[1].CopyNew()
	if out[3], err = a.xor(x[2], x[7]); err != nil {
		return out, fmt.Errorf("xtime out3: %w", err)
	}
	if out[4], err = a.xor(x[3], x[7]); err != nil {
		return out, fmt.Errorf("xtime out4: %w", err)
	}
	out[5] = x[4].CopyNew()
	out[6] = x[5].CopyNew()
	out[7] = x[6].CopyNew()
	return out, nil
}

// AddRoundKey XORs a round key byte by byte
func (a *Evaluator) AddRoundKey(st, rk StateHE) (StateHE, error) {
	var out StateHE
	for b := 0; b < 16; b++ {
		ob, err := a.xorByte(st[b], rk[b])
		if err != nil {
			return out, fmt.Errorf("AddRoundKey byte %d: %w", b, err)
		}
		out[b] = ob
	}
	return out, nil
}

// ShiftRowsHE is pointer shuffling only (free)
func ShiftRowsHE(st StateHE) StateHE {
	var out StateHE
	for r := 0; r < 4; r++ {
		for c := 0; c < 4; c++ {
			out[r+4*c] = st[r+4*((c+r)%4)]
		}
	}
	return out
}

// mixColumnsV1 is the original MixColumns (byte-level xtime + xorBytes), depth 4, with the
// output uniformized by squaring.
//
// Kept for comparison, the active version is MixColumnsV2
// (V2, bit-level tree, depth 3). Per column [b0,b1,b2,b3]:
//
//	D0 = xtime(b0 ^ b1) ^ b1 ^ b2 ^ b3
//	D1 = xtime(b1 ^ b2) ^ b2 ^ b3 ^ b0
//	D2 = xtime(b2 ^ b3) ^ b3 ^ b0 ^ b1
//	D3 = xtime(b3 ^ b0) ^ b0 ^ b1 ^ b2
func (a *Evaluator) mixColumnsV1(st StateHE) (StateHE, error) {
	var out StateHE
	for c := 0; c < 4; c++ {
		i := 4 * c
		b0, b1, b2, b3 := st[i], st[i+1], st[i+2], st[i+3]

		col := [4]struct {
			pair0, pair1 ByteHE
			t0, t1, t2   ByteHE
		}{
			{b0, b1, b1, b2, b3}, // D0
			{b1, b2, b2, b3, b0}, // D1
			{b2, b3, b3, b0, b1}, // D2
			{b3, b0, b0, b1, b2}, // D3
		}

		for d := 0; d < 4; d++ {
			u, err := a.xorByte(col[d].pair0, col[d].pair1)
			if err != nil {
				return out, fmt.Errorf("MixColumns col %d D%d xor pair: %w", c, d, err)
			}
			t, err := a.xtime(u)
			if err != nil {
				return out, fmt.Errorf("MixColumns col %d D%d xtime: %w", c, d, err)
			}
			D, err := a.xorBytes(t, col[d].t0, col[d].t1, col[d].t2)
			if err != nil {
				return out, fmt.Errorf("MixColumns col %d D%d reduce: %w", c, d, err)
			}
			out[i+d] = D
		}
	}

	// Level uniformization after MixColumns, the paths consumes different depths (xtime+reduce
	// = 4 levels, cloned paths = 3), so the state is not "flat".
	cts := make([]*rlwe.Ciphertext, 0, 128)
	for b := 0; b < 16; b++ {
		for i := 0; i < 8; i++ {
			cts = append(cts, out[b][i])
		}
	}
	if a.noTrick {
		// Standard leveled uniformization: DropLevel only, no squaring.
		utils.FlattenLevels(a.eval, cts)
	} else if err := utils.FlattenBitLevels(a.eval, cts); err != nil {
		return out, fmt.Errorf("MixColumns level uniformization: %w", err)
	}

	return out, nil
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
	return reduceBalanced(cts, a.xor)
}

// MixColumnsV2: each output bit is a balanced XOR tree over its fresh input leaves.
// At most 7 leaves => depth ceil(log2 7) = 3
func (a *Evaluator) MixColumnsV2(st StateHE) (StateHE, error) {
	var out StateHE
	for c := 0; c < 4; c++ {
		base := 4 * c
		for d := 0; d < 4; d++ {
			var ob ByteHE
			for j := 0; j < 8; j++ {
				leaves := mcLeavesTable[d][j]
				cts := make([]*rlwe.Ciphertext, len(leaves))
				for n, lf := range leaves {
					cts[n] = st[base+lf.off][lf.bit]
				}
				r, err := a.xorTree(cts)
				if err != nil {
					return out, fmt.Errorf("MixColumnsV2 col %d D%d bit %d: %w", c, d, j, err)
				}
				ob[j] = r
			}
			out[base+d] = ob
		}
	}
	return out, nil
}

// SubBytes applies the selected SubByte version (1..4) to all 16 bytes of the state.
//
// SubBytes is the one round op whose trick/no-trick behaviour is chosen by the version rather
// than by the evaluator's noTrick flag: V1 builds the S-box monomials with the squaring product
// (utils.MulBits) and flattens by squaring (utils.FlattenBitLevels).
func (a *Evaluator) SubBytes(st StateHE, version int) (StateHE, error) {
	var out StateHE
	if a.noTrick && version == 1 {
		return out, fmt.Errorf("SubBytes: version 1 is the squaring-trick S-box (MulBits + FlattenBitLevels); use -subbytes 2, 3 or 4 in no-trick mode")
	}
	var sub func(ByteHE) (ByteHE, error)
	switch version {
	case 1:
		sub = a.SubByteV1
	case 2:
		sub = a.SubByteV2
	case 3:
		sub = a.SubByteV3
	case 4:
		sub = a.SubByteV4
	default:
		return out, fmt.Errorf("SubBytes: unknown version %d (want 1..4)", version)
	}
	for b := 0; b < 16; b++ {
		t0 := time.Now()
		ob, err := sub(st[b])
		if err != nil {
			fmt.Println()
			return out, fmt.Errorf("SubBytes byte %d (v%d): %w", b, version, err)
		}
		fmt.Printf(" [SubBytes v%d] byte %2d/16 done (%s)\n", version, b+1, time.Since(t0).Round(time.Millisecond))
		out[b] = ob
	}
	return out, nil
}

// MixColumns dispatches to the selected MixColumns version
func (a *Evaluator) MixColumns(st StateHE, version int) (StateHE, error) {
	switch version {
	case 1:
		return a.mixColumnsV1(st)
	case 2:
		return a.MixColumnsV2(st)
	default:
		return StateHE{}, fmt.Errorf("MixColumns: unknown version %d (want 1..2)", version)
	}
}
