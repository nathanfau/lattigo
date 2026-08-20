package blt

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"

	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charctx"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charenc"
)

type LTStrategy int

const (
	LTBSGS LTStrategy = iota
	LTNaive
	LTBSGSWithGiantStep
)

type CompileOptions struct {
	Strategy  LTStrategy
	GiantStep int
}

type CompiledTransform struct {
	In, Out       charenc.BlockSpec
	LT            *lintrans.LinearTransformation
	DiagPlaintext *rlwe.Plaintext
	BiasPlaintext *rlwe.Plaintext
	GaloisEls     []uint64
	InputLevel    int
	Raw           RawCompiled
}

// Exactly one of LT or DiagPlaintext is non-nil. The DiagPlaintext path
// sidesteps MultiplyByDiagMatrix never initialising the output when the only
// non-zero LT diagonal is index 0.
type RawCompiled struct {
	LT            *lintrans.LinearTransformation
	DiagPlaintext *rlwe.Plaintext
	BiasPlaintext *rlwe.Plaintext
	GaloisEls     []uint64
	InputLevel    int
	OutSlots      int
	InSlots       int
}

// ltScaleAt is hefloat.Parameters.GetScalingFactor(DefaultScale, DefaultScale, level), which v6
// does not expose: the product of the primes a rescale consumes at that level, so the LT output
// lands back at DefaultScale once rescaled.
func ltScaleAt(params ckks.Parameters, level int) rlwe.Scale {
	b := params.DefaultScale()
	Q := params.Q()
	for i := 0; i < params.LevelsConsumedPerRescaling(); i++ {
		b = b.Mul(rlwe.NewScale(Q[level-i]))
	}
	return b.Div(params.DefaultScale())
}

// Both LT plaintext and bias plaintext are pre-scaled so the post-rescale
// output lands at DefaultScale.
func CompileMatrix(matrix [][]complex128, bias []complex128, outSlots, inSlots int, ctx *charctx.Context, inputLevel int) (*RawCompiled, error) {
	if inputLevel < ctx.Params.LevelsConsumedPerRescaling() {
		return nil, fmt.Errorf("blt.CompileMatrix: inputLevel %d below LevelsConsumedPerRescaling %d", inputLevel, ctx.Params.LevelsConsumedPerRescaling())
	}
	slots := ctx.Params.MaxSlots()
	if inSlots > slots || outSlots > slots {
		return nil, fmt.Errorf("blt.CompileMatrix: block does not fit in %d slots", slots)
	}
	if len(matrix) != outSlots {
		return nil, fmt.Errorf("blt.CompileMatrix: matrix has %d rows, expected %d", len(matrix), outSlots)
	}
	if bias != nil && len(bias) != outSlots {
		return nil, fmt.Errorf("blt.CompileMatrix: bias length %d, expected %d", len(bias), outSlots)
	}

	diagonals := lintrans.Diagonals[complex128]{}
	diagonalOnly := true
	for i := 0; i < outSlots; i++ {
		if len(matrix[i]) != inSlots {
			return nil, fmt.Errorf("blt.CompileMatrix: row %d has %d cols, expected %d", i, len(matrix[i]), inSlots)
		}
		for k := 0; k < inSlots; k++ {
			v := matrix[i][k]
			if v == 0 {
				continue
			}
			if k != i {
				diagonalOnly = false
			}
			d := (k - i) % slots
			if d < 0 {
				d += slots
			}
			row, ok := diagonals[d]
			if !ok {
				row = make([]complex128, slots)
				diagonals[d] = row
			}
			row[i] = v
		}
	}

	return CompileDiagonals(diagonals, diagonalOnly, bias, outSlots, inSlots, ctx, inputLevel)
}

func CompileDiagonals(diagonals lintrans.Diagonals[complex128], diagonalOnly bool, bias []complex128, outSlots, inSlots int, ctx *charctx.Context, inputLevel int) (*RawCompiled, error) {
	return CompileDiagonalsWithOptions(diagonals, diagonalOnly, bias, outSlots, inSlots, ctx, inputLevel, CompileOptions{})
}

// Lattigo convention: each diagonals[d] has ctx.Params.MaxSlots() entries and
// output slot i receives diagonals[d][i] * input slot (i+d).
func CompileDiagonalsWithOptions(diagonals lintrans.Diagonals[complex128], diagonalOnly bool, bias []complex128, outSlots, inSlots int, ctx *charctx.Context, inputLevel int, opts CompileOptions) (*RawCompiled, error) {
	if inputLevel < ctx.Params.LevelsConsumedPerRescaling() {
		return nil, fmt.Errorf("blt.CompileDiagonalsWithOptions: inputLevel %d below LevelsConsumedPerRescaling %d", inputLevel, ctx.Params.LevelsConsumedPerRescaling())
	}
	slots := ctx.Params.MaxSlots()
	if inSlots > slots || outSlots > slots {
		return nil, fmt.Errorf("blt.CompileDiagonalsWithOptions: block does not fit in %d slots", slots)
	}
	if bias != nil && len(bias) != outSlots {
		return nil, fmt.Errorf("blt.CompileDiagonalsWithOptions: bias length %d, expected %d", len(bias), outSlots)
	}
	for d, row := range diagonals {
		if len(row) != slots {
			return nil, fmt.Errorf("blt.CompileDiagonalsWithOptions: diagonal %d has %d slots, expected %d", d, len(row), slots)
		}
	}

	ltScale := ltScaleAt(ctx.Params, inputLevel)

	if diagonalOnly {
		diagValues := make([]complex128, slots)
		if row, ok := diagonals[0]; ok {
			copy(diagValues, row)
		}
		pt := ckks.NewPlaintext(ctx.Params, inputLevel)
		pt.Scale = ltScale
		if err := ctx.Encoder.Encode(diagValues, pt); err != nil {
			return nil, fmt.Errorf("blt.CompileDiagonalsWithOptions: encode diagonal plaintext: %w", err)
		}
		var biasPt *rlwe.Plaintext
		if bias != nil && hasNonZero(bias) {
			outLevel := inputLevel - ctx.Params.LevelsConsumedPerRescaling()
			padded := make([]complex128, slots)
			copy(padded, bias)
			biasPt = ckks.NewPlaintext(ctx.Params, outLevel)
			if err := ctx.Encoder.Encode(padded, biasPt); err != nil {
				return nil, fmt.Errorf("blt.CompileDiagonalsWithOptions: encode bias: %w", err)
			}
		}
		return &RawCompiled{
			DiagPlaintext: pt,
			BiasPlaintext: biasPt,
			InputLevel:    inputLevel,
			OutSlots:      outSlots,
			InSlots:       inSlots,
		}, nil
	}

	giantStep, err := giantStepFor(opts)
	if err != nil {
		return nil, fmt.Errorf("blt.CompileDiagonalsWithOptions: %w", err)
	}
	indexes := diagonals.DiagonalsIndexList()
	if giantStep == 0 && len(indexes) < 2 {
		// Lattigo's BSGS auto giant-step requires at least two non-zero diagonals.
		giantStep = -1
	}

	ltParams := lintrans.Parameters{
		DiagonalsIndexList:        indexes,
		LevelQ:                    inputLevel,
		LevelP:                    ctx.Params.MaxLevelP(),
		Scale:                     ltScale,
		LogDimensions:             ctx.Params.LogMaxDimensions(),
		LogBabyStepGiantStepRatio: giantStep,
	}
	lt := lintrans.NewTransformation(ctx.Params, ltParams)
	if err := lintrans.Encode(ctx.Encoder, diagonals, lt); err != nil {
		return nil, fmt.Errorf("blt.CompileDiagonalsWithOptions: encode LT: %w", err)
	}
	galEls := lintrans.GaloisElements(ctx.Params, ltParams)

	var biasPt *rlwe.Plaintext
	if bias != nil && hasNonZero(bias) {
		outLevel := inputLevel - ctx.Params.LevelsConsumedPerRescaling()
		padded := make([]complex128, slots)
		copy(padded, bias)
		biasPt = ckks.NewPlaintext(ctx.Params, outLevel)
		if err := ctx.Encoder.Encode(padded, biasPt); err != nil {
			return nil, fmt.Errorf("blt.CompileDiagonalsWithOptions: encode bias: %w", err)
		}
	}

	return &RawCompiled{
		LT:            &lt,
		BiasPlaintext: biasPt,
		GaloisEls:     galEls,
		InputLevel:    inputLevel,
		OutSlots:      outSlots,
		InSlots:       inSlots,
	}, nil
}

// giantStepFor returns v6's LogBabyStepGiantStepRatio. LTNaive and LTBSGS carry over unchanged
// (-1 disables BSGS, 0 lets lattigo pick). LTBSGSWithGiantStep does NOT: Pro7ech's GiantStep is
// the giant step n1 itself, whereas v6 takes log2(n1/n2). The value is passed through as that
// log ratio; it is a benchmark-only knob, off by default.
func giantStepFor(opts CompileOptions) (int, error) {
	switch opts.Strategy {
	case LTNaive:
		return -1, nil
	case LTBSGS:
		return 0, nil
	case LTBSGSWithGiantStep:
		if opts.GiantStep <= 0 {
			return 0, fmt.Errorf("LTBSGSWithGiantStep requires GiantStep > 0, got %d", opts.GiantStep)
		}
		return opts.GiantStep, nil
	default:
		return 0, fmt.Errorf("unknown LT strategy %d", opts.Strategy)
	}
}

func Compile(tr Transform, ctx *charctx.Context, inputLevel int) (*CompiledTransform, error) {
	if tr.Out.Slots == 0 || tr.In.Slots == 0 {
		return nil, fmt.Errorf("blt.Compile: degenerate spec In.Slots=%d Out.Slots=%d", tr.In.Slots, tr.Out.Slots)
	}
	raw, err := CompileMatrix(tr.Matrix, tr.Bias, tr.Out.Slots, tr.In.Slots, ctx, inputLevel)
	if err != nil {
		return nil, err
	}
	compiled := &CompiledTransform{
		In:            tr.In,
		Out:           tr.Out,
		LT:            raw.LT,
		DiagPlaintext: raw.DiagPlaintext,
		BiasPlaintext: raw.BiasPlaintext,
		GaloisEls:     raw.GaloisEls,
		InputLevel:    raw.InputLevel,
		Raw:           *raw,
	}
	return compiled, nil
}

func hasNonZero(v []complex128) bool {
	for _, z := range v {
		if z != 0 {
			return true
		}
	}
	return false
}
