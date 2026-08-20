package lut

import (
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/blt"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charctx"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charenc"
)

// Only the unary path is ported. Upstream this file also carries the tensor evaluator
// (EvalBinary, EvalArity3, EvalArity4, EvalTensor) and the two helpers that fan its stages out
// over the blt worker pool, parallelApplyRaw and parallelMulRelin. They all take a
// *CompiledTensor, which lives in the unported multivariate.go, so they come as a block if the
// multivariate LUTs are ever needed.

type Evaluator struct {
	Ctx *charctx.Context
	BLT *blt.Evaluator
}

func NewEvaluator(ctx *charctx.Context) *Evaluator {
	return &Evaluator{Ctx: ctx, BLT: blt.NewEvaluator(ctx)}
}

func NewEvaluatorWithWorkerCapacity(ctx *charctx.Context, capacity int) *Evaluator {
	return &Evaluator{Ctx: ctx, BLT: blt.NewEvaluatorWithCapacity(ctx, capacity)}
}

func (e *Evaluator) EvalUnary(x charenc.CipherBlock, lut *CompiledUnary) (charenc.CipherBlock, error) {
	return e.BLT.Apply(x, lut.Transform)
}
