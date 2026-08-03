package aes

import "github.com/tuneinsight/lattigo/v6/core/rlwe"

var (
	shiftRowsSrc  = [8]int{0, 5, 2, 7, 4, 1, 6, 3}
	shiftRowsSwap = [8]bool{false, false, true, true, false, true, true, false}
)

// ShiftRows applies AES ShiftRows at the Algo1 pause, where the state is already split into its
// Re / Im streams: pure pointer moves, no homomorphic operation at all.
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
