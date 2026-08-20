// Package whaes builds the AES round on the packed Walsh-Hadamard byte layout, where SubBytes,
// ShiftRows and the GF(2^8) multiplications of MixColumns fuse into one block linear transform
// per output column,
//
//	M_r' = sum_alpha (P_{r',alpha} . P_pi) tensor A_alpha
//
// It is the WH counterpart of nathanfau/aes: same cipher, but the S-box is a matrix instead of a
// monomial circuit, because the character basis linearises any map of a byte. The cleartext
// reference and the GF(2^8) arithmetic come from aes, only the homomorphic circuit differs.
//
// One byte is one WH block of t = 2^8 slots, so a state takes 16*256 = 2^12 slots, in the STANDARD
// ring for now. This file holds A_alpha, the S-box matrix.
package whaes

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/blt"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charenc"
)

// T is the alphabet size of one AES byte, t = 2^8.
const T = 256

// Spec is the block spec of one AES byte. not reduced
func Spec() charenc.BlockSpec {
	return charenc.BlockSpec{
		Alphabet: charenc.Alphabet{Modulus: T, Kind: charenc.WH},
		Reduced:  false,
	}
}

// Codec is the WH codec of one AES byte.
func Codec() (charenc.Codec, error) {
	return charenc.NewCodec(Spec())
}

// SBoxMap is S_r : m -> r . S(m). A non-zero r makes it bijective, hence the zero bias.
func SBoxMap(r byte) func(int) int {
	return func(m int) int { return int(aes.GFMul(r, aes.SBox(byte(m)))) }
}

// SBoxTransform is the 256x256 matrix A_r with WH(r . S(m)) = A_r . WH(m) for every m, obtained by
// interpolating S_r on the WH basis.
func SBoxTransform(r byte) (blt.Transform, error) {
	if r == 0 {
		return blt.Transform{}, fmt.Errorf("whaes.SBoxTransform: r must be non-zero")
	}
	codec, err := Codec()
	if err != nil {
		return blt.Transform{}, fmt.Errorf("whaes.SBoxTransform: codec: %w", err)
	}
	tr, err := blt.CompileUnary(codec, codec, SBoxMap(r))
	if err != nil {
		return blt.Transform{}, fmt.Errorf("whaes.SBoxTransform: interpolate: %w", err)
	}
	return tr, nil
}

// SBoxMatrix is SBoxTransform's matrix alone.
func SBoxMatrix(r byte) ([][]complex128, error) {
	tr, err := SBoxTransform(r)
	if err != nil {
		return nil, err
	}
	return tr.Matrix, nil
}
