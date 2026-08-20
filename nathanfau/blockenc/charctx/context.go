package charctx

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"

	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charenc"
)

// Gen is bumped whenever Evaluator is rebuilt via WithKey; shallow-copied
// workers must be discarded when it advances.
type Context struct {
	Params      ckks.Parameters
	Encoder     *ckks.Encoder
	Encryptor   *rlwe.Encryptor
	Decryptor   *rlwe.Decryptor
	KeyGen      *rlwe.KeyGenerator
	SecretKey   *rlwe.SecretKey
	EvK         *rlwe.MemEvaluationKeySet
	Evaluator   *ckks.Evaluator
	LTEvaluator *lintrans.Evaluator
	Gen         int

	// GaloisLevelQ caps the level the Galois keys are generated at.
	GaloisLevelQ *int
}

func NewContext(pl ckks.ParametersLiteral) (*Context, error) {
	params, err := ckks.NewParametersFromLiteral(pl)
	if err != nil {
		return nil, fmt.Errorf("charctx: invalid parameters: %w", err)
	}
	ctx := &Context{Params: params}
	ctx.KeyGen = ckks.NewKeyGenerator(params)
	ctx.SecretKey = ctx.KeyGen.GenSecretKeyNew()
	rlk := ctx.KeyGen.GenRelinearizationKeyNew(ctx.SecretKey)
	ctx.EvK = rlwe.NewMemEvaluationKeySet(rlk)
	ctx.Encoder = ckks.NewEncoder(params)
	ctx.Encryptor = rlwe.NewEncryptor(params, ctx.SecretKey)
	ctx.Decryptor = rlwe.NewDecryptor(params, ctx.SecretKey)
	ctx.Evaluator = ckks.NewEvaluator(params, ctx.EvK)
	ctx.LTEvaluator = lintrans.NewEvaluator(ctx.Evaluator)
	return ctx, nil
}

// EnsureGaloisKeys registers any missing Galois keys and rebuilds the
// Evaluator/LTEvaluator via WithKey so Lattigo's automorphism-index cache
// covers the new elements.
func (ctx *Context) EnsureGaloisKeys(elements []uint64) {
	anyMissing := false
	for _, e := range elements {
		if _, present := ctx.EvK.GaloisKeys[e]; !present {
			anyMissing = true
			break
		}
	}
	if !anyMissing {
		return
	}
	seen := make(map[uint64]struct{}, len(elements))
	missing := make([]uint64, 0, len(elements))
	for _, e := range elements {
		if _, dup := seen[e]; dup {
			continue
		}
		seen[e] = struct{}{}
		if _, present := ctx.EvK.GaloisKeys[e]; !present {
			missing = append(missing, e)
		}
	}
	if len(missing) == 0 {
		return
	}
	var evkParams []rlwe.EvaluationKeyParameters
	if ctx.GaloisLevelQ != nil {
		evkParams = append(evkParams, rlwe.EvaluationKeyParameters{LevelQ: ctx.GaloisLevelQ})
	}
	keys := ctx.KeyGen.GenGaloisKeysNew(missing, ctx.SecretKey, evkParams...)
	for _, k := range keys {
		ctx.EvK.GaloisKeys[k.GaloisElement] = k
	}
	ctx.Evaluator = ctx.Evaluator.WithKey(ctx.EvK)
	ctx.LTEvaluator = lintrans.NewEvaluator(ctx.Evaluator)
	ctx.Gen++
}

func (ctx *Context) EncodeBlock(block charenc.PlainBlock, level int) (*rlwe.Plaintext, error) {
	slots := ctx.Params.MaxSlots()
	if len(block.Values) > slots {
		return nil, fmt.Errorf("charctx: block has %d slots, ciphertext supports %d", len(block.Values), slots)
	}
	padded := make([]complex128, slots)
	copy(padded, block.Values)
	pt := ckks.NewPlaintext(ctx.Params, level)
	if err := ctx.Encoder.Encode(padded, pt); err != nil {
		return nil, fmt.Errorf("charctx: encode: %w", err)
	}
	return pt, nil
}

func (ctx *Context) DecodeBlock(pt *rlwe.Plaintext, spec charenc.BlockSpec, layout charenc.Layout) (charenc.PlainBlock, error) {
	slots := ctx.Params.MaxSlots()
	values := make([]complex128, slots)
	if err := ctx.Encoder.Decode(pt, values); err != nil {
		return charenc.PlainBlock{}, fmt.Errorf("charctx: decode: %w", err)
	}
	return charenc.PlainBlock{
		Spec:   spec,
		Layout: layout,
		Values: values[:spec.Slots],
	}, nil
}

func (ctx *Context) EncryptBlock(block charenc.PlainBlock, level int) (charenc.CipherBlock, error) {
	pt, err := ctx.EncodeBlock(block, level)
	if err != nil {
		return charenc.CipherBlock{}, err
	}
	ct := ckks.NewCiphertext(ctx.Params, 1, level)
	if err := ctx.Encryptor.Encrypt(pt, ct); err != nil {
		return charenc.CipherBlock{}, fmt.Errorf("charctx: encrypt: %w", err)
	}
	layout := block.Layout
	if layout.BlockSlots == 0 {
		layout = charenc.Layout{Kind: charenc.BlockInOneCiphertext, BlockSlots: block.Spec.Slots, Stride: 1}
	}
	return charenc.CipherBlock{Spec: block.Spec, Layout: layout, CTs: []*rlwe.Ciphertext{ct}}, nil
}

func (ctx *Context) DecryptBlock(cb charenc.CipherBlock) (charenc.PlainBlock, error) {
	if len(cb.CTs) != 1 {
		return charenc.PlainBlock{}, fmt.Errorf("charctx: expected one ciphertext, got %d", len(cb.CTs))
	}
	pt := ctx.Decryptor.DecryptNew(cb.CTs[0])
	return ctx.DecodeBlock(pt, cb.Spec, cb.Layout)
}

// DefaultParameters is a testing-only CKKS set; not secure for production.
func DefaultParameters() ckks.ParametersLiteral {
	return ckks.ParametersLiteral{
		LogN:            10,
		LogQ:            []int{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{60},
		LogDefaultScale: 40,
	}
}
