package cleaning

import (
	"fmt"
	"math"
	"math/rand"
	"testing"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func newCleaningContext(t *testing.T) (ckks.Parameters, *ckks.Encoder, *rlwe.Encryptor, *rlwe.Decryptor, *ckks.Evaluator) {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{55, 40, 40, 40, 40},
		LogP:            []int{55},
		LogDefaultScale: 40,
	})
	if err != nil {
		t.Fatalf("params: %v", err)
	}
	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	evk := rlwe.NewMemEvaluationKeySet(rlk)

	ecd := ckks.NewEncoder(params)
	enc := rlwe.NewEncryptor(params, pk)
	dec := rlwe.NewDecryptor(params, sk)
	eval := ckks.NewEvaluator(params, evk)

	debug.EncStd, debug.DecStd, debug.ParamsStd = ecd, dec, params

	return params, ecd, enc, dec, eval
}

// TestCleaning runs Cleaning on noisy bits (an exact 0/1 plus a perturbation) at several
// noise levels, and checks it computes p(x)=3x^2-2x^3 exactly, is scale-neutral, and
// consumes 2 levels. It logs input vs output error to show the contraction (out ~ 3*in^2).
func TestCleaning(t *testing.T) {
	params, ecd, enc, dec, eval := newCleaningContext(t)
	n := params.MaxSlots()

	for _, noise := range []float64{1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9} {
		t.Run(fmt.Sprintf("noise=%.0e", noise), func(t *testing.T) {
			rng := rand.New(rand.NewSource(1))
			bit := make([]float64, n)
			in := make([]float64, n)
			for s := 0; s < n; s++ {
				b := float64(rng.Intn(2))
				bit[s] = b
				in[s] = b + (rng.Float64()*2-1)*noise
			}

			pt := ckks.NewPlaintext(params, params.MaxLevel())
			if err := ecd.Encode(in, pt); err != nil {
				t.Fatalf("encode: %v", err)
			}
			ct, err := enc.EncryptNew(pt)
			if err != nil {
				t.Fatalf("encrypt: %v", err)
			}
			W := ct.Scale
			inLevel := ct.Level()

			want := make([]float64, n)
			for s := 0; s < n; s++ {
				x := in[s]
				want[s] = 3*x*x - 2*x*x*x
			}

			debug.DbgSlotStd(fmt.Sprintf("noisy bits (%.0e) =", noise), ct)
			debug.DbgChain("Chain before cleaning :", eval, ct)
			debug.PrintPrecStd("prec in  :", want, ct)

			out, err := Cleaning(params, eval, ct)
			if err != nil {
				t.Fatalf("Cleaning: %v", err)
			}
			debug.DbgSlotStd("cleaned bits =", out)
			debug.DbgChain("Chain after cleaning :", eval, out)
			debug.PrintPrecStd("prec out :", want, out)

			got := make([]complex128, out.Slots())
			if err := ecd.Decode(dec.DecryptNew(out), got); err != nil {
				t.Fatalf("decode: %v", err)
			}

			maxErr := 0.0
			var inErr, outErr float64
			for s := 0; s < n; s++ {
				if e := math.Abs(real(got[s]) - want[s]); e > maxErr {
					maxErr = e
				}
				inErr += math.Abs(in[s] - bit[s])
				outErr += math.Abs(real(got[s]) - bit[s])
			}
			inErr /= float64(n)
			outErr /= float64(n)
			t.Logf("noise=%.0e: value err vs p=%.3e | mean err to bit: in=%.3e -> out=%.3e", noise, maxErr, inErr, outErr)

			if maxErr > 1e-4 {
				t.Fatalf("Cleaning value error %.3e too large", maxErr)
			}
			if rel := math.Abs(out.Scale.Float64()-W.Float64()) / W.Float64(); rel > 1e-9 {
				t.Fatalf("Cleaning not scale-neutral: out 2^%.4f, want 2^%.4f", math.Log2(out.Scale.Float64()), math.Log2(W.Float64()))
			}
			if d := inLevel - out.Level(); d != 2 {
				t.Fatalf("Cleaning consumed %d levels, want 2", d)
			}
		})
	}
}
