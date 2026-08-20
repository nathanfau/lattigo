package blt

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charenc"
)

// Affine block linear transform A * omega + b mapping In coordinates to Out.
type Transform struct {
	In, Out charenc.BlockSpec
	Matrix  [][]complex128 // Out.Slots x In.Slots
	Bias    []complex128   // Out.Slots
	Layout  charenc.Layout
}

func CompileUnary(in, out charenc.Codec, f func(int) int) (Transform, error) {
	t := in.Spec().Alphabet.Modulus
	if t <= 0 {
		return Transform{}, fmt.Errorf("blt: empty input alphabet")
	}
	values := make([][]complex128, t)
	for m := 0; m < t; m++ {
		values[m] = out.EncodeValue(f(m))
	}
	A, b, err := in.Interpolate(values, out.Spec())
	if err != nil {
		return Transform{}, err
	}
	return Transform{
		In:     in.Spec(),
		Out:    out.Spec(),
		Matrix: A,
		Bias:   b,
	}, nil
}

func Apply(in charenc.PlainBlock, T Transform) (charenc.PlainBlock, error) {
	if in.Spec != T.In {
		return charenc.PlainBlock{}, fmt.Errorf("blt.Apply: input spec %+v does not match transform In %+v", in.Spec, T.In)
	}
	if len(in.Values) != T.In.Slots {
		return charenc.PlainBlock{}, fmt.Errorf("blt.Apply: input has %d slots, expected %d", len(in.Values), T.In.Slots)
	}
	if len(T.Bias) != T.Out.Slots {
		return charenc.PlainBlock{}, fmt.Errorf("blt.Apply: bias has %d entries, expected %d", len(T.Bias), T.Out.Slots)
	}
	out := make([]complex128, T.Out.Slots)
	for j := 0; j < T.Out.Slots; j++ {
		if len(T.Matrix[j]) != T.In.Slots {
			return charenc.PlainBlock{}, fmt.Errorf("blt.Apply: matrix row %d has %d cols, expected %d", j, len(T.Matrix[j]), T.In.Slots)
		}
		sum := T.Bias[j]
		for k := 0; k < T.In.Slots; k++ {
			sum += T.Matrix[j][k] * in.Values[k]
		}
		out[j] = sum
	}
	return charenc.PlainBlock{
		Spec:   T.Out,
		Layout: in.Layout,
		Values: out,
	}, nil
}
