package whaes

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/nathanfau/blockpack"
)

// The 13 small 16x16 byte patterns of the round: P_pi and the twelve P_{r',alpha}. The state is
// numbered FIPS-197, byte 4c+r

const Bytes = 16

// MixAlpha is the MixColumns matrix, alpha[r][r'].
var MixAlpha = [4][4]byte{
	{2, 3, 1, 1},
	{1, 2, 3, 1},
	{1, 1, 2, 3},
	{3, 1, 1, 2},
}

// Alphas are the three GF(2^8) coefficients MixColumns needs, hence the three big matrices.
var Alphas = [3]byte{1, 2, 3}

// ShiftRowsPerm returns perm with perm[i] the input byte feeding output byte i, i.e. pi^{-1}.
// It mirrors aes.ShiftRows.
func ShiftRowsPerm() (perm [Bytes]int) {
	for r := 0; r < 4; r++ {
		for c := 0; c < 4; c++ {
			perm[r+4*c] = r + 4*((c+r)%4)
		}
	}
	return perm
}

// ShiftRowsMatrix is P_pi: P[i][j] = 1 iff j = pi^{-1}(i).
func ShiftRowsMatrix() [][]complex128 {
	perm := ShiftRowsPerm()
	P := Zeros(Bytes, Bytes)
	for i, j := range perm {
		P[i][j] = 1
	}
	return P
}

// MixPartition is P_{r',alpha} = {r : alpha_{r,r'} = alpha}. For each r' the three sets partition
// {0,1,2,3}, which is what makes the sum of the three alpha terms well defined slot by slot.
func MixPartition(rPrime int, alpha byte) []int {
	var out []int
	for r := 0; r < 4; r++ {
		if MixAlpha[r][rPrime] == alpha {
			out = append(out, r)
		}
	}
	return out
}

// MixSelect is the 16x16 matrix P_{r',alpha}: byte 4c+r' goes to every 4c+r with r in the
// partition, everything else is zeroed.
func MixSelect(rPrime int, alpha byte) ([][]complex128, error) {
	if rPrime < 0 || rPrime > 3 {
		return nil, fmt.Errorf("whaes.MixSelect: r' = %d, want 0..3", rPrime)
	}
	P := Zeros(Bytes, Bytes)
	for _, r := range MixPartition(rPrime, alpha) {
		for c := 0; c < 4; c++ {
			P[4*c+r][4*c+rPrime] = 1
		}
	}
	return P, nil
}

// FusedPattern is Q_{r',alpha} = P_{r',alpha} . P_pi, so that
//
//	M_r' = sum_alpha Q_{r',alpha} tensor A_alpha
//
// carries SubBytes, ShiftRows, the GF(2^8) multiplication and the MixColumns dispatch at once.
func FusedPattern(rPrime int, alpha byte) ([][]complex128, error) {
	sel, err := MixSelect(rPrime, alpha)
	if err != nil {
		return nil, err
	}
	return Mul(sel, ShiftRowsMatrix())
}

// KronSum is sum_i patterns[i] tensor blocks[i], kept factored: M_r' densely would be a 4096x4096
// complex matrix, 268 MB, and every use here reads straight off the factors.
type KronSum struct {
	Patterns [][][]complex128 // 16x16 each
	Blocks   [][][]complex128 // 256x256 each
}

// RoundKronSum is M_r' in factored form.
func RoundKronSum(rPrime int, sbox [3][][]complex128) (KronSum, error) {
	ks := KronSum{}
	for i, alpha := range Alphas {
		q, err := FusedPattern(rPrime, alpha)
		if err != nil {
			return KronSum{}, err
		}
		ks.Patterns = append(ks.Patterns, q)
		ks.Blocks = append(ks.Blocks, sbox[i])
	}
	return ks, nil
}

// Slots is blockpack.WHSlotVec for one state, in complex128 because that is what the round
// matrices work in. Only the element type changes here.
func Slots(slots int, s [Bytes]byte) []complex128 {
	vals := blockpack.WHSlotVec(slots, [][16]byte{s})
	out := make([]complex128, len(vals))
	for i, v := range vals {
		out[i] = complex(v, 0)
	}
	return out
}

// Apply evaluates the sum on a vector of Bytes*T slots, block by block.
func (ks KronSum) Apply(x []complex128) ([]complex128, error) {
	if len(x) != Bytes*T {
		return nil, fmt.Errorf("whaes.KronSum.Apply: %d slots, want %d", len(x), Bytes*T)
	}
	out := make([]complex128, Bytes*T)
	for k, pat := range ks.Patterns {
		A := ks.Blocks[k]
		for i := 0; i < Bytes; i++ {
			for j := 0; j < Bytes; j++ {
				w := pat[i][j]
				if w == 0 {
					continue
				}
				in := x[j*T : (j+1)*T]
				dst := out[i*T : (i+1)*T]
				for u := 0; u < T; u++ {
					acc := complex128(0)
					row := A[u]
					for v := 0; v < T; v++ {
						acc += row[v] * in[v]
					}
					dst[u] += w * acc
				}
			}
		}
	}
	return out, nil
}

// DiagonalValues returns the assembled matrix in lattigo's diagonal form, diag[d][i] sending input
// slot (i+d) mod slots to output slot i. The second return says whether diagonal 0 carries it all,
// which lets blt take its single-plaintext shortcut.
func (ks KronSum) DiagonalValues(slots int) (map[int][]complex128, bool) {
	diag := map[int][]complex128{}
	diagonalOnly := true
	for k, pat := range ks.Patterns {
		A := ks.Blocks[k]
		for ib := 0; ib < len(pat); ib++ {
			for jb := 0; jb < len(pat[ib]); jb++ {
				w := pat[ib][jb]
				if w == 0 {
					continue
				}
				for u := 0; u < len(A); u++ {
					i := ib*T + u
					for v := 0; v < len(A[u]); v++ {
						if A[u][v] == 0 {
							continue
						}
						j := jb*T + v
						if i != j {
							diagonalOnly = false
						}
						d := (j - i) % slots
						if d < 0 {
							d += slots
						}
						row, ok := diag[d]
						if !ok {
							row = make([]complex128, slots)
							diag[d] = row
						}
						row[i] += w * A[u][v]
					}
				}
			}
		}
	}
	return diag, diagonalOnly
}

// Diagonals is the set of non-zero diagonal indices, what drives the homomorphic cost of a BLT.
func (ks KronSum) Diagonals(slots int) map[int]struct{} {
	diag := map[int]struct{}{}
	for k, pat := range ks.Patterns {
		A := ks.Blocks[k]

		off := map[int]struct{}{}
		for u := 0; u < len(A); u++ {
			for v := 0; v < len(A[u]); v++ {
				if A[u][v] != 0 {
					off[v-u] = struct{}{}
				}
			}
		}
		// The pattern is 16x16 for a full state, 1x1 when a single block is measured on its own.
		shift := map[int]struct{}{}
		for i := 0; i < len(pat); i++ {
			for j := 0; j < len(pat[i]); j++ {
				if pat[i][j] != 0 {
					shift[j-i] = struct{}{}
				}
			}
		}
		for s := range shift {
			for e := range off {
				d := (T*s + e) % slots
				if d < 0 {
					d += slots
				}
				diag[d] = struct{}{}
			}
		}
	}
	return diag
}
