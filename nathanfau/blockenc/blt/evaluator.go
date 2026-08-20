package blt

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"

	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charctx"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charenc"
)

// A ckks.Evaluator at logN=15 holds roughly 30 MB of RNS buffers, so this also bounds memory.
const DefaultWorkerPoolCapacity = 8

// Workers carry the ctx.Gen at which they were created; stale ones are
// discarded when the context's evaluator is rebuilt by EnsureGaloisKeys.
type Evaluator struct {
	Ctx     *charctx.Context
	workers chan *Worker
	slots   chan struct{}
	cap     int
}

// Workers are not safe for concurrent use; each goroutine must acquire its
// own via GetWorker.
type Worker struct {
	Eval   *ckks.Evaluator
	LTEval *lintrans.Evaluator
	gen    int
}

func NewEvaluator(ctx *charctx.Context) *Evaluator {
	return NewEvaluatorWithCapacity(ctx, DefaultWorkerPoolCapacity)
}

func NewEvaluatorWithCapacity(ctx *charctx.Context, capacity int) *Evaluator {
	if capacity < 1 {
		capacity = 1
	}
	return &Evaluator{
		Ctx:     ctx,
		workers: make(chan *Worker, capacity),
		slots:   make(chan struct{}, capacity),
		cap:     capacity,
	}
}

func (e *Evaluator) Capacity() int {
	return e.cap
}

// newWorker builds a private evaluator rather than shallow-copying the context's: v6 dropped
// ShallowCopy, and WithKey shares its temporary buffers, so it cannot be used concurrently.
func (e *Evaluator) newWorker() *Worker {
	sub := ckks.NewEvaluator(e.Ctx.Params, e.Ctx.EvK)
	return &Worker{
		Eval:   sub,
		LTEval: lintrans.NewEvaluator(sub),
		gen:    e.Ctx.Gen,
	}
}

func (e *Evaluator) GetWorker() *Worker {
	for {
		select {
		case w := <-e.workers:
			if w.gen != e.Ctx.Gen {
				e.releaseWorker()
				continue
			}
			return w
		default:
		}

		select {
		case e.slots <- struct{}{}:
			return e.newWorker()
		case w := <-e.workers:
			if w.gen != e.Ctx.Gen {
				e.releaseWorker()
				continue
			}
			return w
		}
	}
}

func (e *Evaluator) PutWorker(w *Worker) {
	if w.gen != e.Ctx.Gen {
		e.releaseWorker()
		return
	}
	select {
	case e.workers <- w:
	default:
		e.releaseWorker()
	}
}

func (e *Evaluator) releaseWorker() {
	select {
	case <-e.slots:
	default:
	}
}

func (e *Evaluator) Apply(in charenc.CipherBlock, ct *CompiledTransform) (charenc.CipherBlock, error) {
	if in.Spec != ct.In {
		return charenc.CipherBlock{}, fmt.Errorf("blt.Apply: input spec %+v does not match compiled In %+v", in.Spec, ct.In)
	}
	if len(in.CTs) != 1 {
		return charenc.CipherBlock{}, fmt.Errorf("blt.Apply: expected 1 ciphertext for block layout, got %d", len(in.CTs))
	}
	e.Ctx.EnsureGaloisKeys(ct.GaloisEls)
	w := e.GetWorker()
	defer e.PutWorker(w)
	out, err := e.applyRawWith(w, in.CTs[0], &ct.Raw)
	if err != nil {
		return charenc.CipherBlock{}, err
	}
	return charenc.CipherBlock{
		Spec:   ct.Out,
		Layout: in.Layout,
		CTs:    []*rlwe.Ciphertext{out},
	}, nil
}

func (e *Evaluator) ApplyRaw(in *rlwe.Ciphertext, ct *RawCompiled) (*rlwe.Ciphertext, error) {
	e.Ctx.EnsureGaloisKeys(ct.GaloisEls)
	w := e.GetWorker()
	defer e.PutWorker(w)
	return e.applyRawWith(w, in, ct)
}

// Caller must pre-ensure the required Galois keys.
func (e *Evaluator) ApplyRawWith(w *Worker, in *rlwe.Ciphertext, ct *RawCompiled) (*rlwe.Ciphertext, error) {
	return e.applyRawWith(w, in, ct)
}

func (e *Evaluator) applyRawWith(w *Worker, in *rlwe.Ciphertext, ct *RawCompiled) (*rlwe.Ciphertext, error) {
	if in.Level() < ct.InputLevel {
		return nil, fmt.Errorf("blt.applyRawWith: input level %d below compiled InputLevel %d", in.Level(), ct.InputLevel)
	}
	out := ckks.NewCiphertext(e.Ctx.Params, 1, ct.InputLevel)

	if ct.DiagPlaintext != nil {
		if err := w.Eval.Mul(in, ct.DiagPlaintext, out); err != nil {
			return nil, fmt.Errorf("blt.applyRawWith: diag mul: %w", err)
		}
	} else {
		// The BSGS/naive split and the hoisting buffers are internal to v6's evaluator: it
		// decomposes the input and pre-rotates on its own, driven by the LT's own N1. The
		// upstream version drove that by hand to cache the pre-rotations across calls; that
		// cache is dropped here, nothing else changes.
		if err := w.LTEval.Evaluate(in, *ct.LT, out); err != nil {
			return nil, fmt.Errorf("blt.applyRawWith: LT evaluate: %w", err)
		}
	}
	if err := w.Eval.Rescale(out, out); err != nil {
		return nil, fmt.Errorf("blt.applyRawWith: rescale: %w", err)
	}
	if ct.BiasPlaintext != nil {
		if err := w.Eval.Add(out, ct.BiasPlaintext, out); err != nil {
			return nil, fmt.Errorf("blt.applyRawWith: bias add: %w", err)
		}
	}
	return out, nil
}
