package charenc

import "github.com/tuneinsight/lattigo/v6/core/rlwe"

type LayoutKind int

const (
	BlockInOneCiphertext LayoutKind = iota
	CoordinatePerCiphertext
)

type Layout struct {
	Kind       LayoutKind
	BlockSlots int
	Stride     int
	Offset     int
}

type PlainBlock struct {
	Spec   BlockSpec
	Layout Layout
	Values []complex128
}

// CipherBlock: CTs has length 1 for BlockInOneCiphertext, Spec.Slots for
// CoordinatePerCiphertext.
type CipherBlock struct {
	Spec   BlockSpec
	Layout Layout
	CTs    []*rlwe.Ciphertext
}
