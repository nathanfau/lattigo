package lut

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/blt"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charctx"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charenc"
)

type CompiledUnary struct {
	Table     UnaryTable
	Transform *blt.CompiledTransform
}

func CompileUnary(table UnaryTable, ctx *charctx.Context, inputLevel int) (*CompiledUnary, error) {
	in, err := charenc.NewCodec(table.In)
	if err != nil {
		return nil, fmt.Errorf("lut.CompileUnary: input codec: %w", err)
	}
	out, err := charenc.NewCodec(table.Out)
	if err != nil {
		return nil, fmt.Errorf("lut.CompileUnary: output codec: %w", err)
	}
	tr, err := blt.CompileUnary(in, out, table.Eval)
	if err != nil {
		return nil, fmt.Errorf("lut.CompileUnary: BLT: %w", err)
	}
	cct, err := blt.Compile(tr, ctx, inputLevel)
	if err != nil {
		return nil, fmt.Errorf("lut.CompileUnary: BLT compile: %w", err)
	}
	ctx.EnsureGaloisKeys(cct.GaloisEls)
	return &CompiledUnary{Table: table, Transform: cct}, nil
}
