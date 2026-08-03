// target case: 2^15 blocks, logN = 16
// no problem when logN < 16 in tests, it adapts
//
// AES state (row-major, Byte = 4*row + col + 1):
//
//	[ Byte1  ][ Byte2  ][ Byte3  ][ Byte4  ]
//	[ Byte5  ][ Byte6  ][ Byte7  ][ Byte8  ]
//	[ Byte9  ][ Byte10 ][ Byte11 ][ Byte12 ]
//	[ Byte13 ][ Byte14 ][ Byte15 ][ Byte16 ]
//
// 64 ciphertexts indexed [group][bit], 8 bits per group (low half byte | high half byte):
//
//	ct  1..8   (group 0) : [ Byte1  | Byte3  ]
//	ct  9..16  (group 1) : [ Byte5  | Byte7  ]
//	ct 17..24  (group 2) : [ Byte9  | Byte11 ]
//	ct 25..32  (group 3) : [ Byte13 | Byte15 ]
//	ct 33..40  (group 4) : [ Byte2  | Byte4  ]
//	ct 41..48  (group 5) : [ Byte6  | Byte8  ]
//	ct 49..56  (group 6) : [ Byte10 | Byte12 ]
//	ct 57..64  (group 7) : [ Byte14 | Byte16 ]

package blockpack

import (
	"fmt"
	"math/rand"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// Packed is a block-batched, bit-sliced AES state: 64 ciphertexts indexed [group][bit].
type Packed [8][8]*rlwe.Ciphertext

// Group g's low half carries AES byte lowByte[g], its high half byte highByte[g]. This pairing is
// the layout the HE round functions keep invariant.
var (
	lowByte  = [8]int{0, 4, 8, 12, 1, 5, 9, 13}
	highByte = [8]int{2, 6, 10, 14, 3, 7, 11, 15}
)

// Capacity is how many AES blocks fit in one Packed state: MaxSlots/2, each block taking one slot
// per half.
func Capacity(params ckks.Parameters) int {
	return params.MaxSlots() / 2
}

// Cts flattens a Packed into its 64 ciphertexts, in [group][bit] order.
func (p Packed) Cts() []*rlwe.Ciphertext {
	cts := make([]*rlwe.Ciphertext, 0, 64)
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			cts = append(cts, p[g][b])
		}
	}
	return cts
}

// SlotVec builds the cleartext slot vector of ciphertext [g][b]: block s on slot s (low half) and
// slot half+s (high half), unused slots left at zero. It is what Encrypt encodes.
func SlotVec(params ckks.Parameters, blocks [][16]byte, g, b int) []float64 {
	slots := params.MaxSlots()
	half := slots / 2
	vals := make([]float64, slots)
	for s, blk := range blocks {
		vals[s] = float64((blk[lowByte[g]] >> uint(b)) & 1)
		vals[half+s] = float64((blk[highByte[g]] >> uint(b)) & 1)
	}
	return vals
}

// RandomBlocks returns Capacity(params) independent random blocks. Distinct blocks are what
// exercise the packing: a replicated block passes even if the layout mixes slots up.
func RandomBlocks(params ckks.Parameters, rng *rand.Rand) [][16]byte {
	blocks := make([][16]byte, Capacity(params))
	for i := range blocks {
		for b := 0; b < 16; b++ {
			blocks[i][b] = byte(rng.Intn(256))
		}
	}
	return blocks
}

// Check decrypts p and compares it to want block by block, naming the first mismatch.
func Check(params ckks.Parameters, ecd *ckks.Encoder, dec *rlwe.Decryptor, p Packed, want [][16]byte) error {
	got, worst, err := Decrypt(params, ecd, dec, p)
	if err != nil {
		return fmt.Errorf("decrypt: %w", err)
	}
	if len(got) != len(want) {
		return fmt.Errorf("decrypted %d blocks, want %d", len(got), len(want))
	}
	bad, first := 0, -1
	for s := range want {
		if got[s] != want[s] {
			bad++
			if first < 0 {
				first = s
			}
		}
	}
	if bad != 0 {
		return fmt.Errorf("%d/%d blocks wrong (worst bit err %.3f); block %d got = %x, want = %x",
			bad, len(want), worst, first, got[first], want[first])
	}
	return nil
}

// WantSlots returns the expected slot vectors of the 64 ciphertexts, in the order Cts uses: what
// utils.BitDistance expects alongside it.
func WantSlots(params ckks.Parameters, blocks [][16]byte) [][]float64 {
	out := make([][]float64, 0, 64)
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			out = append(out, SlotVec(params, blocks, g, b))
		}
	}
	return out
}

// Encrypt packs up to Capacity(params) blocks into a Packed state at the given level.
func Encrypt(params ckks.Parameters, ecd *ckks.Encoder, enc *rlwe.Encryptor, blocks [][16]byte, level int) (Packed, error) {
	var p Packed
	slots := params.MaxSlots()
	half := slots / 2
	if len(blocks) > half {
		return p, fmt.Errorf("blockpack.Encrypt: %d blocks exceed capacity %d (= MaxSlots/2)", len(blocks), half)
	}

	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			vals := SlotVec(params, blocks, g, b)
			pt := ckks.NewPlaintext(params, level)
			if err := ecd.Encode(vals, pt); err != nil {
				return p, fmt.Errorf("blockpack.Encrypt encode [g=%d][b=%d]: %w", g, b, err)
			}
			ct := ckks.NewCiphertext(params, 1, level)
			if err := enc.Encrypt(pt, ct); err != nil {
				return p, fmt.Errorf("blockpack.Encrypt encrypt [g=%d][b=%d]: %w", g, b, err)
			}
			p[g][b] = ct
		}
	}
	return p, nil
}

// Decrypt reverses Encrypt, rounding each decrypted value to the nearest bit. It also returns the
// worst |value - bit| seen over every slot: a small margin means clean bits, a margin near 0.5
// means the blocks below were read out of noise and are meaningless.
func Decrypt(params ckks.Parameters, ecd *ckks.Encoder, dec *rlwe.Decryptor, p Packed) ([][16]byte, float64, error) {
	slots := params.MaxSlots()
	half := slots / 2
	out := make([][16]byte, half)
	buf := make([]complex128, slots)
	worst := 0.0

	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			if p[g][b] == nil {
				return nil, 0, fmt.Errorf("blockpack.Decrypt: nil ciphertext [g=%d][b=%d]", g, b)
			}
			if err := ecd.Decode(dec.DecryptNew(p[g][b]), buf); err != nil {
				return nil, 0, fmt.Errorf("blockpack.Decrypt decode [g=%d][b=%d]: %w", g, b, err)
			}
			if d := utils.BitDistanceOf(buf, nil, 0).Worst; d > worst {
				worst = d
			}
			for s := 0; s < half; s++ {
				if real(buf[s]) >= 0.5 {
					out[s][lowByte[g]] |= 1 << uint(b)
				}
				if real(buf[half+s]) >= 0.5 {
					out[s][highByte[g]] |= 1 << uint(b)
				}
			}
		}
	}
	return out, worst, nil
}
