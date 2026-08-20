package lut

import "github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charenc"

type UnaryTable struct {
	In   charenc.BlockSpec
	Out  charenc.BlockSpec
	Eval func(x int) int
}

type MultiTable struct {
	In   []charenc.BlockSpec
	Out  charenc.BlockSpec
	Eval func(xs []int) int
}
