// Walsh-Hadamard packing of an AES state, the layout of nathanfau/whtransciphering.
//
// One byte m is one WH block of t = 2^8 slots, WH(m) = ((-1)^{s.m})_{s in F_2^8}, so a 16-byte
// state takes 16*256 = 4096 slots and sits in ONE ciphertext. No reordering: byte i of the block
// goes to WH block i, unlike the bit-sliced packing above which pairs bytes into groups.
package blockpack

import (
	"fmt"
	"math"
	"math/bits"
	"math/rand"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const (
	// WHBytes is the number of bytes of an AES state.
	WHBytes = 16
	// WHT is the alphabet of one byte, t = 2^8, hence the slots of one WH block.
	WHT = 256
	// WHSlotsPerBlock is what one packed state costs: 16 bytes of 256 characters.
	WHSlotsPerBlock = WHBytes * WHT
)

// WHCapacity is how many AES blocks fit in slots.
func WHCapacity(slots int) int {
	return slots / WHSlotsPerBlock
}

// WHSlot is the slot carrying character s of byte g for block i.
func WHSlot(slots, g, s, i int) int {
	return (g*WHT+s)*WHCapacity(slots) + i
}

// WHSlotVec builds the cleartext slot vector of a packed state.
func WHSlotVec(slots int, blocks [][16]byte) []float64 {
	capacity := WHCapacity(slots)
	vals := make([]float64, slots)
	for i := range vals {
		vals[i] = 1
	}
	for i, blk := range blocks {
		for g := 0; g < WHBytes; g++ {
			m := blk[g]
			for s := 0; s < WHT; s++ {
				if bits.OnesCount8(uint8(s)&m)&1 == 1 {
					vals[(g*WHT+s)*capacity+i] = -1
				}
			}
		}
	}
	return vals
}

// WHReadSlots reads the blocks back from a decoded slot vector, and returns the worst
// | |value| - 1 | over every slot. A clean slot is +-1, so that distance says how far the encoding
// has drifted
//
// The decoder only reads the sign of the 8 characters s = 2^j of each byte because the rest of the
// block is redundant, but the drift is measured over all of them.
func WHReadSlots(buf []complex128) ([][16]byte, float64) {
	capacity := WHCapacity(len(buf))

	worst := 0.0
	for k := 0; k < capacity*WHSlotsPerBlock; k++ {
		if d := math.Abs(math.Abs(real(buf[k])) - 1); d > worst {
			worst = d
		}
	}

	out := make([][16]byte, capacity)
	for i := 0; i < capacity; i++ {
		for g := 0; g < WHBytes; g++ {
			var m byte
			for j := 0; j < 8; j++ {
				if real(buf[(g*WHT+(1<<uint(j)))*capacity+i]) < 0 {
					m |= 1 << uint(j)
				}
			}
			out[i][g] = m
		}
	}
	return out, worst
}

// WHRandomBlocks returns WHCapacity(slots) independent random blocks. Distinct blocks are what
// exercise the interleaving: a replicated block passes even if the stride is wrong.
func WHRandomBlocks(slots int, rng *rand.Rand) [][16]byte {
	blocks := make([][16]byte, WHCapacity(slots))
	for i := range blocks {
		for b := 0; b < WHBytes; b++ {
			blocks[i][b] = byte(rng.Intn(256))
		}
	}
	return blocks
}

// WHEncrypt packs up to WHCapacity blocks into one ciphertext at the given level.
func WHEncrypt(params ckks.Parameters, ecd *ckks.Encoder, enc *rlwe.Encryptor, blocks [][16]byte, level int) (*rlwe.Ciphertext, error) {
	slots := params.MaxSlots()
	capacity := WHCapacity(slots)
	if capacity == 0 {
		return nil, fmt.Errorf("blockpack.WHEncrypt: %d slots, one WH state needs %d", slots, WHSlotsPerBlock)
	}
	if len(blocks) > capacity {
		return nil, fmt.Errorf("blockpack.WHEncrypt: %d blocks exceed capacity %d", len(blocks), capacity)
	}

	pt := ckks.NewPlaintext(params, level)
	if err := ecd.Encode(WHSlotVec(slots, blocks), pt); err != nil {
		return nil, fmt.Errorf("blockpack.WHEncrypt encode: %w", err)
	}
	ct := ckks.NewCiphertext(params, 1, level)
	if err := enc.Encrypt(pt, ct); err != nil {
		return nil, fmt.Errorf("blockpack.WHEncrypt encrypt: %w", err)
	}
	return ct, nil
}

// WHDecrypt reverses WHEncrypt: decrypt, decode, then WHReadSlots.
func WHDecrypt(params ckks.Parameters, ecd *ckks.Encoder, dec *rlwe.Decryptor, ct *rlwe.Ciphertext) ([][16]byte, float64, error) {
	if ct == nil {
		return nil, 0, fmt.Errorf("blockpack.WHDecrypt: nil ciphertext")
	}
	slots := params.MaxSlots()
	if WHCapacity(slots) == 0 {
		return nil, 0, fmt.Errorf("blockpack.WHDecrypt: %d slots, one WH state needs %d", slots, WHSlotsPerBlock)
	}

	buf := make([]complex128, slots)
	if err := ecd.Decode(dec.DecryptNew(ct), buf); err != nil {
		return nil, 0, fmt.Errorf("blockpack.WHDecrypt decode: %w", err)
	}
	blocks, worst := WHReadSlots(buf)
	return blocks, worst, nil
}

// WHCheck decrypts ct and compares it to want block by block, naming the first mismatch.
func WHCheck(params ckks.Parameters, ecd *ckks.Encoder, dec *rlwe.Decryptor, ct *rlwe.Ciphertext, want [][16]byte) error {
	got, worst, err := WHDecrypt(params, ecd, dec, ct)
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
		return fmt.Errorf("%d/%d blocks wrong (worst |value|-1 err %.3f); block %d got = %x, want = %x",
			bad, len(want), worst, first, got[first], want[first])
	}
	return nil
}
