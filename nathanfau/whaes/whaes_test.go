package whaes

//	go test ./nathanfau/whaes/ -v

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"sort"
	"testing"

	"github.com/tuneinsight/lattigo/v6/circuits/common/lintrans"
	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockpack"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
)

// sboxMatrices is A_1, A_2, A_3, which both tests need.
func sboxMatrices(t *testing.T) [3][][]complex128 {
	t.Helper()
	var out [3][][]complex128
	for i, alpha := range Alphas {
		m, err := SBoxMatrix(alpha)
		if err != nil {
			t.Fatalf("SBoxMatrix(%d): %v", alpha, err)
		}
		out[i] = m
	}
	return out
}

// TestFusedRound : the four M_r' applied to the encoded state, XORed together slot-wise,
// must give SubBytes then ShiftRows then MixColumns.
func TestFusedRound(t *testing.T) {

	sboxMats := sboxMatrices(t)
	rng := rand.New(rand.NewSource(1))

	fmt.Printf("================ fused round, %d states ================\n", 8)

	for iter := 0; iter < 8; iter++ {
		var st [Bytes]byte
		for i := range st {
			st[i] = byte(rng.Intn(256))
		}

		want := st
		aes.SubBytes(want[:])
		aes.ShiftRows(want[:])
		aes.MixColumns(want[:])

		// T_r' = M_r' . ct, then XOR the four, i.e. multiply slot-wise.
		x := Slots(Bytes*T, st)
		acc := make([]complex128, Bytes*T)
		for i := range acc {
			acc[i] = 1
		}
		for rPrime := 0; rPrime < 4; rPrime++ {
			ks, err := RoundKronSum(rPrime, sboxMats)
			if err != nil {
				t.Fatalf("RoundKronSum(%d): %v", rPrime, err)
			}
			Tr, err := ks.Apply(x)
			if err != nil {
				t.Fatalf("Apply(%d): %v", rPrime, err)
			}
			for i := range acc {
				acc[i] *= Tr[i]
			}
		}

		blocks, _ := blockpack.WHReadSlots(acc)
		if got := blocks[0]; got != want {
			t.Fatalf("iter %d: got %x, want %x (state %x)", iter, got, want, st)
		}
	}
	fmt.Printf("  8/8 states matchs SubBytes+ShiftRows+MixColumns\n")
}

// diagStats is one matrix measured: what it costs to store it, to apply it, and how much of what
// is stored is actually non-zero.
type diagStats struct {
	diags   int     // non-zero diagonals: one encoded plaintext each
	rots    int     // key-switches under BSGS, baby steps plus giant steps
	nonzero int     // non-zero coefficients over all diagonals
	minAbs  float64 // smallest and largest |coefficient|
	maxAbs  float64
}

// fill is the share of the stored slots that carry something. A diagonal costs a full plaintext
// whatever its density, so what this measures is waste.
func (s diagStats) fill() float64 { return float64(s.nonzero) / float64(s.diags*Bytes*T) }

// measure assembles the diagonals, reads the coefficients off them, then runs lattigo's own BSGS
// split to get the rotation count the evaluator would pay.
func measure(ks KronSum, slots int) diagStats {
	vals, _ := ks.DiagonalValues(slots)
	st := diagStats{diags: len(vals), minAbs: math.Inf(1)}

	idx := make([]int, 0, len(vals))
	for d, row := range vals {
		idx = append(idx, d)
		for _, v := range row {
			if v == 0 {
				continue
			}
			st.nonzero++
			a := cmplx.Abs(v)
			if a < st.minAbs {
				st.minAbs = a
			}
			if a > st.maxAbs {
				st.maxAbs = a
			}
		}
	}
	sort.Ints(idx)

	n1 := lintrans.FindBestBSGSRatio(idx, slots, 1)
	index, _, _ := lintrans.BSGSIndex(idx, slots, n1)
	baby, giant := map[int]struct{}{}, 0
	for j, is := range index {
		if j != 0 {
			giant++
		}
		for _, i := range is {
			if i != 0 {
				baby[i] = struct{}{}
			}
		}
	}
	st.rots = len(baby) + giant

	vals = nil
	return st
}

// TestRoundDiagonals measures what the round costs homomorphically
func TestRoundDiagonals(t *testing.T) {

	sboxMats := sboxMatrices(t)
	slots := Bytes * T

	single := KronSum{Patterns: [][][]complex128{Identity(1)}, Blocks: [][][]complex128{sboxMats[0]}}
	subBytes := KronSum{Patterns: [][][]complex128{Identity(Bytes)}, Blocks: [][][]complex128{sboxMats[0]}}
	fusedSR := KronSum{Patterns: [][][]complex128{ShiftRowsMatrix()}, Blocks: [][][]complex128{sboxMats[0]}}

	leads := []string{"A_1 (one byte)", "I16 x A_1 (SubBytes)", "Ppi x A_1 (+ShiftRows)"}
	all := []diagStats{measure(single, slots), measure(subBytes, slots), measure(fusedSR, slots)}

	for rPrime := 0; rPrime < 4; rPrime++ {
		ks, err := RoundKronSum(rPrime, sboxMats)
		if err != nil {
			t.Fatalf("RoundKronSum(%d): %v", rPrime, err)
		}
		leads = append(leads, fmt.Sprintf("M_%d (+MixColumns)", rPrime))
		all = append(all, measure(ks, slots))
	}

	rows := make([][]string, len(all))
	for i, s := range all {
		rows[i] = []string{
			fmt.Sprintf("%d", s.diags),
			fmt.Sprintf("%d", s.rots),
			fmt.Sprintf("%d", s.nonzero),
			fmt.Sprintf("%.1f%%", 100*s.fill()),
			fmt.Sprintf("%.4f", s.minAbs),
			fmt.Sprintf("%.4f", s.maxAbs),
		}
	}

	fmt.Printf("the round matrices over %d slots\n%s", slots,
		utils.BoxTable("matrix",
			[]utils.TableGroup{
				{Name: "cost", Cols: []string{"diagonals", "rotations"}},
				{Name: "coefficients", Cols: []string{"non-zero", "fill", "min |v|", "max |v|"}},
			},
			leads, rows))

	// Tensoring by I16 is free: the block shifts stay at {0}, so the whole state costs what one
	// byte costs. Fusing ShiftRows is not: P_pi brings four distinct block shifts.
	if all[1].diags != all[0].diags {
		t.Errorf("I16 x A_1 has %d diagonals, A_1 alone has %d: tensoring by the identity should be free", all[1].diags, all[0].diags)
	}
	if all[2].diags <= all[1].diags {
		t.Errorf("Ppi x A_1 has %d diagonals, I16 x A_1 has %d: fusing ShiftRows should cost", all[2].diags, all[1].diags)
	}
}
