// The homomorphic side of the package: the evaluator every round function hangs off, and the
// selectors that tie the circuits together. One file per round function next to this one --
// sbox_he.go, shiftrows_he.go, mixcolumns_he.go, xor_he.go -- and aes.go for the cleartext AES the
// tests check against.
package aes

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
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

// XorFunc is a two-ciphertext XOR circuit: XorSq and XorNoSq both fit.
type XorFunc func(x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error)

// XorPlainFunc is the same against a cleartext bit plane: XorSqPlain and XorNoSqPlain both fit.
type XorPlainFunc func(x *rlwe.Ciphertext, bits []float64) (*rlwe.Ciphertext, error)

// XorKind names a circuit, so one value picks both flavours at once.
type XorKind int

const (
	NoSq XorKind = iota // x+y-2xy, the pipeline default and the zero value
	Sq                  // (x-y)^2
)

func (k XorKind) String() string {
	if k == Sq {
		return "sq ((x-y)^2)"
	}
	return "nosq (x+y-2xy)"
}

// ParseXorKind reads the name a flag carries.
func ParseXorKind(s string) (XorKind, error) {
	switch s {
	case "sq", "square":
		return Sq, nil
	case "nosq", "arith":
		return NoSq, nil
	}
	return 0, fmt.Errorf("unknown XOR %q, want \"sq\" or \"nosq\"", s)
}

// Of and PlainOf return the XOR circuits of the given kind.
func (a *Evaluator) Of(k XorKind) XorFunc {
	if k == Sq {
		return a.XorSq
	}
	return a.XorNoSq
}

func (a *Evaluator) PlainOf(k XorKind) XorPlainFunc {
	if k == Sq {
		return a.XorSqPlain
	}
	return a.XorNoSqPlain
}

// SubByte applies the S-box variant named by version (1..3), cheapest circuit last.
func (a *Evaluator) SubByte(inp ByteHE, version int) (ByteHE, error) {
	switch version {
	case 1:
		return a.SubByteV1(inp)
	case 2:
		return a.SubByteV2(inp)
	case 3:
		return a.SubByteV3(inp)
	default:
		return ByteHE{}, fmt.Errorf("SubByte: unknown version %d (want 1..3)", version)
	}
}
