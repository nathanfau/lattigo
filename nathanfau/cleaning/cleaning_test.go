package cleaning

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// expSizes are the error magnitudes tested by the benchmark, every slot gets |err| ~ 2^-exp.
var expSizes = []int{2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30}

type cleanFn struct {
	name string
	fn   func(ckks.Parameters, *ckks.Evaluator, *rlwe.Ciphertext) (*rlwe.Ciphertext, error)
}

var cleanFns = []cleanFn{
	{"Cleaning", Cleaning},
	{"Smoother", SmootherCleaning},
	{"VerySmoother", VerySmootherCleaning},
}

type cell struct {
	precBits  float64
	worstBits float64
	dur       time.Duration
}

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

	return params, ecd, enc, dec, eval
}

// noisyBits returns n bits perturbed by an error of the order of 2^-exp
func noisyBits(rng *rand.Rand, n, exp int) (bit, in []float64) {
	noise := math.Exp2(-float64(exp))
	jitter := jitterFrac(exp)
	bit = make([]float64, n)
	in = make([]float64, n)
	for s := 0; s < n; s++ {
		b := float64(rng.Intn(2))
		bit[s] = b
		sign := float64(rng.Intn(2)*2 - 1)
		mag := noise * (1 + (rng.Float64()*2-1)*jitter)
		in[s] = b + sign*mag
	}
	return bit, in
}

func jitterFrac(exp int) float64 {
	j := float64(exp) / 32
	if j > 0.9 {
		j = 0.9
	}
	return j
}

// TestCleaning sweeps every cleaning polynomial over a range of input error sizes
// and prints a single table comparing the output precision to the bit and the run time.
func TestCleaning(t *testing.T) {
	params, ecd, enc, dec, eval := newCleaningContext(t)
	n := params.MaxSlots()

	res := make([][]cell, len(expSizes))
	for i := range res {
		res[i] = make([]cell, len(cleanFns))
	}
	// Starting precision to the bit, per error size (average over the slots, and worst slot).
	inPrec := make([]float64, len(expSizes))
	inWorst := make([]float64, len(expSizes))

	for ei, exp := range expSizes {
		rng := rand.New(rand.NewSource(int64(exp)))
		bit, in := noisyBits(rng, n, exp)

		pt := ckks.NewPlaintext(params, params.MaxLevel())
		if err := ecd.Encode(in, pt); err != nil {
			t.Fatalf("encode: %v", err)
		}
		ct, err := enc.EncryptNew(pt)
		if err != nil {
			t.Fatalf("encrypt: %v", err)
		}
		st, err := utils.BitDistanceCt(ecd, dec, ct, bit, 0)
		if err != nil {
			t.Fatalf("input precision (exp=2^-%d): %v", exp, err)
		}
		inPrec[ei], inWorst[ei] = st.AvgPrec, -math.Log2(st.Worst)

		for fi, cf := range cleanFns {
			t0 := time.Now()
			out, err := cf.fn(params, eval, ct)
			dur := time.Since(t0)
			if err != nil {
				t.Fatalf("%s (exp=2^-%d): %v", cf.name, exp, err)
			}
			st, err := utils.BitDistanceCt(ecd, dec, out, bit, 0)
			if err != nil {
				t.Fatalf("%s (exp=2^-%d): %v", cf.name, exp, err)
			}
			res[ei][fi] = cell{st.AvgPrec, -math.Log2(st.Worst), dur}
		}
	}

	printCleaningTable(res, inPrec, inWorst)

	// Sanity: in the clearly-contracting regime, cleaning must not lose precision to the bit.
	for ei, exp := range expSizes {
		if exp > 12 {
			continue
		}
		for fi, cf := range cleanFns {
			if out := res[ei][fi].precBits; out < inPrec[ei]-0.5 {
				t.Errorf("%s regressed at exp=2^-%d: in=%.2f b -> out=%.2f b", cf.name, exp, inPrec[ei], out)
			}
		}
	}
}

// printCleaningTable renders the benchmark as a box-drawing grid (via utils.BoxTable): one row per
// error size.
// Each precision (avg prec and worst prec) is followed by its gain over the input, and the run time.
func printCleaningTable(res [][]cell, inPrec, inWorst []float64) {
	groups := make([]utils.TableGroup, len(cleanFns))
	for i, cf := range cleanFns {
		groups[i] = utils.TableGroup{Name: cf.name, Cols: []string{"avg", "diff avg", "worst", "diff worst", "time"}}
	}

	leads := make([]string, len(expSizes))
	rows := make([][]string, len(expSizes))
	for ei := range expSizes {
		leads[ei] = fmt.Sprintf("%.2f/%.2f", inPrec[ei], inWorst[ei])
		row := make([]string, 0, 5*len(cleanFns))
		for fi := range cleanFns {
			c := res[ei][fi]
			row = append(row,
				fmt.Sprintf("%.2f", c.precBits), fmt.Sprintf("%+.2f", c.precBits-inPrec[ei]),
				fmt.Sprintf("%.2f", c.worstBits), fmt.Sprintf("%+.2f", c.worstBits-inWorst[ei]),
				c.dur.Round(time.Microsecond).String())
		}
		rows[ei] = row
	}

	fmt.Printf("Cleaning table\n%s",
		utils.BoxTable("in avg/worst", groups, leads, rows))
}
