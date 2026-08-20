package cleaning

//	go test ./nathanfau/cleaning/ -v -timeout 0
//	go test ./nathanfau/cleaning/ -run '^TestCleaning$' -v -timeout 0
//	go test ./nathanfau/cleaning/ -run '^TestXorClean$' -v -timeout 0

import (
	"flag"
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
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

// newCleaningContext is the  context of this file
func newCleaningContext(t *testing.T) (ckks.Parameters, *ckks.Encoder, *rlwe.Encryptor, *rlwe.Decryptor, *ckks.Evaluator) {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{50, 38, 38, 38, 38, 38, 38},
		LogP:            []int{50},
		LogDefaultScale: 38,
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

// signExpSizes stops at 22: past that the output of h sits under the noise floor of the 2^38 scale
// and the quadratic law can no longer be read.
var signExpSizes = []int{1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22}

// noisySigns is noisyBits for the {-1,+1} alphabet: a random sign perturbed by |err| ~ 2^-exp,
// with the same jitter, so the slots do not all sit at the same distance and the worst one means
// something.
func noisySigns(rng *rand.Rand, n, exp int) (sign, in []float64) {
	noise := math.Exp2(-float64(exp))
	jitter := jitterFrac(exp)
	sign = make([]float64, n)
	in = make([]float64, n)
	for s := 0; s < n; s++ {
		b := float64(rng.Intn(2)*2 - 1)
		sign[s] = b
		dir := float64(rng.Intn(2)*2 - 1)
		mag := noise * (1 + (rng.Float64()*2-1)*jitter)
		in[s] = b + dir*mag
	}
	return sign, in
}

// TestSignCleaning sweeps SignCleaning over a range of input error sizes and checks the law it
// rests on: h(±1+e) = ±1 - (3e² ± e³)/2, so |e| -> 1.5|e|² and the worst slot goes from p to
// 2p - log2(1.5) bits. One bit better per level than the {0,1} cleaning 3x²-2x³, which loses 1.58.
func TestSignCleaning(t *testing.T) {
	params, ecd, enc, dec, eval := newCleaningContext(t)
	n := params.MaxSlots()

	const lossBits = 0.5849625007211562 // log2(1.5)

	leads := make([]string, len(signExpSizes))
	rows := make([][]string, len(signExpSizes))

	for ei, exp := range signExpSizes {
		rng := rand.New(rand.NewSource(int64(exp)))
		sign, in := noisySigns(rng, n, exp)

		ct := encryptVals(t, params, ecd, enc, in, params.MaxLevel())
		stIn, err := utils.BitDistanceCt(ecd, dec, ct, sign, 0)
		if err != nil {
			t.Fatalf("input (exp=2^-%d): %v", exp, err)
		}

		t0 := time.Now()
		out, err := SignCleaning(params, eval, ct)
		dur := time.Since(t0)
		if err != nil {
			t.Fatalf("SignCleaning (exp=2^-%d): %v", exp, err)
		}
		if used := params.MaxLevel() - out.Level(); used != 2 {
			t.Errorf("SignCleaning consumed %d levels, want 2", used)
		}

		stOut, err := utils.BitDistanceCt(ecd, dec, out, sign, 0)
		if err != nil {
			t.Fatalf("output (exp=2^-%d): %v", exp, err)
		}

		inWorst, outWorst := -math.Log2(stIn.Worst), -math.Log2(stOut.Worst)
		predicted := 2*inWorst - lossBits

		leads[ei] = fmt.Sprintf("%.2f/%.2f", stIn.AvgPrec, inWorst)
		rows[ei] = []string{
			fmt.Sprintf("%.2f", stOut.AvgPrec), fmt.Sprintf("%+.2f", stOut.AvgPrec-stIn.AvgPrec),
			fmt.Sprintf("%.2f", outWorst), fmt.Sprintf("%+.2f", outWorst-inWorst),
			dur.Round(time.Microsecond).String(),
		}

		// The law only holds while the result stays above the noise floor of the scale.
		if exp <= 12 {
			if outWorst < inWorst {
				t.Errorf("SignCleaning regressed at exp=2^-%d: worst %.2f b -> %.2f b", exp, inWorst, outWorst)
			}
			if math.Abs(outWorst-predicted) > 1 {
				t.Errorf("exp=2^-%d: worst slot %.2f b, the quadratic law predicts %.2f b", exp, outWorst, predicted)
			}
		}
	}

	fmt.Printf("SignCleaning, h(x) = (3x-x^3)/2 on +-1\n%s",
		utils.BoxTable("in avg/worst",
			[]utils.TableGroup{{
				Name: "SignCleaning (2 lv)",
				Cols: []string{"avg", "diff avg", "worst", "diff worst", "time"},
			}},
			leads, rows))
}

type xorCircuit = func(ckks.Parameters, *ckks.Evaluator, *rlwe.Ciphertext, *rlwe.Ciphertext) (*rlwe.Ciphertext, error)

// xorFn is one two-operand circuit: name identifies it in failures, role is its column header.
type xorFn struct {
	name   string
	role   string
	levels int
	fn     xorCircuit
}

// xorPick names one of the package aes XOR circuits, so a table can be built once per circuit.
type xorPick func(*aes.Evaluator) aes.XorFunc

var (
	sq   xorPick = func(a *aes.Evaluator) aes.XorFunc { return a.XorSq }
	noSq xorPick = func(a *aes.Evaluator) aes.XorFunc { return a.XorNoSq }
)

func raw(pick xorPick) xorCircuit {
	return func(_ ckks.Parameters, eval *ckks.Evaluator, x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
		return pick(aes.NewEvaluator(eval))(x, y)
	}
}

func after(pick xorPick, clean CleanFunc) xorCircuit {
	return func(params ckks.Parameters, eval *ckks.Evaluator, x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
		return XorThenClean(params, eval, clean, pick(aes.NewEvaluator(eval)), x, y)
	}
}

func beforeBoth(pick xorPick, clean CleanFunc) xorCircuit {
	return func(params ckks.Parameters, eval *ckks.Evaluator, x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
		return CleanThenXor(params, eval, clean, pick(aes.NewEvaluator(eval)), x, y)
	}
}

func beforeOne(pick xorPick, clean CleanFunc) xorCircuit {
	return func(params ckks.Parameters, eval *ckks.Evaluator, x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
		return CleanOneThenXor(params, eval, clean, pick(aes.NewEvaluator(eval)), x, y)
	}
}

// Which polynomial the XOR tables are built around. 'Cleaning' by default
var xorCleanFlag = flag.String("clean", "cleaning",
	`cleaning polynomial for the XOR tables: "cleaning", "smoother", "verysmoother" or "all"`)

func xorCleans(t *testing.T) []Kind {
	t.Helper()
	if *xorCleanFlag == "all" {
		return []Kind{Basic, Smoother, VerySmoother}
	}
	k, err := ParseKind(*xorCleanFlag)
	if err != nil {
		t.Fatalf("-clean: %v", err)
	}
	return []Kind{k}
}

// xorRow is one input error size of the shared sweep
type xorRow struct {
	want    []float64
	cx, cy  *rlwe.Ciphertext         // both noisy, at MaxLevel
	cyClean *rlwe.Ciphertext         // exact bits, at MaxLevel
	cyAt    map[int]*rlwe.Ciphertext // exact bits, at MaxLevel-depth
	inAvg   float64
	inWorst float64
}

// xorFixture is what the four XOR tests share: one keygen and one set of ciphertexts, built on
// first use and reused, so the tables differ by the circuit alone
type xorFixture struct {
	params ckks.Parameters
	ecd    *ckks.Encoder
	dec    *rlwe.Decryptor
	eval   *ckks.Evaluator
	rows   []xorRow
}

var sharedXorFixture *xorFixture

func xorFix(t *testing.T) *xorFixture {
	t.Helper()
	if sharedXorFixture != nil {
		return sharedXorFixture
	}
	params, ecd, enc, dec, eval := newCleaningContext(t)
	n := params.MaxSlots()

	f := &xorFixture{params: params, ecd: ecd, dec: dec, eval: eval, rows: make([]xorRow, len(expSizes))}
	for ei, e := range expSizes {
		rng := rand.New(rand.NewSource(int64(e)))
		bx, nx := noisyBits(rng, n, e)
		by, ny := noisyBits(rng, n, e)

		want := make([]float64, n)
		for s := range want {
			want[s] = float64(int(bx[s]) ^ int(by[s]))
		}
		cx := encryptVals(t, params, ecd, enc, nx, params.MaxLevel())
		cy := encryptVals(t, params, ecd, enc, ny, params.MaxLevel())

		cyClean := encryptVals(t, params, ecd, enc, by, params.MaxLevel())
		cyAt := map[int]*rlwe.Ciphertext{}
		for _, d := range []int{2, 3} {
			cyAt[d] = encryptVals(t, params, ecd, enc, by, params.MaxLevel()-d)
		}

		sx, err := utils.BitDistanceCt(ecd, dec, cx, bx, 0)
		if err != nil {
			t.Fatalf("input x (e=%d): %v", e, err)
		}
		sy, err := utils.BitDistanceCt(ecd, dec, cy, by, 0)
		if err != nil {
			t.Fatalf("input y (e=%d): %v", e, err)
		}
		f.rows[ei] = xorRow{
			want: want, cx: cx, cy: cy, cyClean: cyClean, cyAt: cyAt,
			inAvg:   math.Min(sx.AvgPrec, sy.AvgPrec),
			inWorst: -math.Log2(math.Max(sx.Worst, sy.Worst)),
		}
	}
	sharedXorFixture = f
	return f
}

// sweepCircuits runs circuits over every error size of the shared fixtur
func sweepCircuits(t *testing.T, f *xorFixture, circuits []xorFn,
	yOf func(ci, ei int) *rlwe.Ciphertext) [][]cell {
	t.Helper()

	res := make([][]cell, len(expSizes))
	for ei := range expSizes {
		r := f.rows[ei]
		res[ei] = make([]cell, len(circuits))
		for ci, c := range circuits {
			t0 := time.Now()
			out, err := c.fn(f.params, f.eval, r.cx, yOf(ci, ei))
			dur := time.Since(t0)
			if err != nil {
				t.Fatalf("%s (e=%d): %v", c.name, expSizes[ei], err)
			}
			if got := f.params.MaxLevel() - out.Level(); got != c.levels {
				t.Errorf("%s consumed %d levels, want %d", c.name, got, c.levels)
			}
			st, err := utils.BitDistanceCt(f.ecd, f.dec, out, r.want, 0)
			if err != nil {
				t.Fatalf("%s (e=%d): %v", c.name, expSizes[ei], err)
			}
			res[ei][ci] = cell{st.AvgPrec, -math.Log2(st.Worst), dur}
		}
	}
	return res
}

// runXorTables prints one table per cleaning polynomial
func runXorTables(t *testing.T, xorName string, pick xorPick, cleanOne bool) {
	t.Helper()
	f := xorFix(t)

	for _, k := range xorCleans(t) {
		clean, depth := k.Func(), k.Depth()

		var circuits []xorFn
		var yOf func(ci, ei int) *rlwe.Ciphertext
		title := fmt.Sprintf("%s, %s", k, xorName)

		if cleanOne {
			circuits = []xorFn{
				{xorName + " clean after", "clean after", 1 + depth, after(pick, clean)},
				{xorName + " clean before", "clean before", depth + 1, beforeBoth(pick, clean)},
				{xorName + " clean x only", "clean x only", depth + 1, beforeOne(pick, clean)},
			}
			// Only the route that never cleans y can take it at h(x)'s level
			yOf = func(ci, ei int) *rlwe.Ciphertext {
				if ci == 2 {
					return f.rows[ei].cyAt[depth]
				}
				return f.rows[ei].cyClean
			}
			title += "   (y = fresh round key, exact bits)"
		} else {
			circuits = []xorFn{
				{xorName, "Xor", 1, raw(pick)},
				{xorName + " clean after", "clean after", 1 + depth, after(pick, clean)},
				{xorName + " clean before", "clean before", depth + 1, beforeBoth(pick, clean)},
			}
			yOf = func(_, ei int) *rlwe.Ciphertext { return f.rows[ei].cy }
		}

		printXorTable(title, circuits, sweepCircuits(t, f, circuits, yOf), f.rows)
	}
}

// printXorTable renders one polynomial the way TestCleaning renders its own sweep: one column group
// per circuit, one row per input error size, each precision followed by its gain over the input.
func printXorTable(name string, circuits []xorFn, res [][]cell, rows []xorRow) {
	groups := make([]utils.TableGroup, len(circuits))
	for i, c := range circuits {
		groups[i] = utils.TableGroup{
			Name: fmt.Sprintf("%s (%d lv)", c.role, c.levels),
			Cols: []string{"avg", "diff avg", "worst", "diff worst", "time"},
		}
	}

	leads := make([]string, len(expSizes))
	out := make([][]string, len(expSizes))
	for ei := range expSizes {
		r := rows[ei]
		leads[ei] = fmt.Sprintf("%.2f/%.2f", r.inAvg, r.inWorst)
		row := make([]string, 0, 5*len(circuits))
		for ci := range circuits {
			c := res[ei][ci]
			row = append(row,
				fmt.Sprintf("%.2f", c.precBits), fmt.Sprintf("%+.2f", c.precBits-r.inAvg),
				fmt.Sprintf("%.2f", c.worstBits), fmt.Sprintf("%+.2f", c.worstBits-r.inWorst),
				c.dur.Round(time.Millisecond).String())
		}
		out[ei] = row
	}

	fmt.Printf("XOR cleaned with %s\n%s", name,
		utils.BoxTable("in avg/worst", groups, leads, out))
}

// TestXorCleanByHand asks whether evaluating the composed polynomial by hand beats going through
// the package Cleaning, which routes h via lattigo's polynomial evaluator and its own scale
// bookkeeping. Two tables per XOR circuit: CleanThenXor on two noisy operands, then
// CleanOneThenXor on a y that is a fresh round key.
func TestXorCleanByHand(t *testing.T) {
	const depth = 2 // the by-hand circuits only implement Cleaning
	f := xorFix(t)

	for _, v := range []struct {
		name      string
		pick      xorPick
		both, one xorCircuit
	}{
		{"XorSq", sq, cleanThenXorSqByHand, cleanOneThenXorSqByHand},
		{"XorNoSq", noSq, cleanThenXorNoSqByHand, cleanOneThenXorNoSqByHand},
	} {
		both := []xorFn{
			{v.name + " both composed", "CleanThenXor", depth + 1, beforeBoth(v.pick, Cleaning)},
			{v.name + " both by hand", "by hand", depth + 1, v.both},
		}
		res := sweepCircuits(t, f, both, func(_, ei int) *rlwe.Ciphertext { return f.rows[ei].cy })
		printXorTable("Cleaning, "+v.name+"   CleanThenXor", both, res, f.rows)

		one := []xorFn{
			{v.name + " one composed", "CleanOneThenXor", depth + 1, beforeOne(v.pick, Cleaning)},
			{v.name + " one by hand", "by hand", depth + 1, v.one},
		}
		res = sweepCircuits(t, f, one, func(_, ei int) *rlwe.Ciphertext { return f.rows[ei].cyAt[depth] })
		printXorTable("Cleaning, "+v.name+"   CleanOneThenXor (y = round key)", one, res, f.rows)
	}
}

func TestXorCleanSq(t *testing.T)      { runXorTables(t, "XorSq", sq, false) }
func TestXorCleanNoSq(t *testing.T)    { runXorTables(t, "XorNoSq", noSq, false) }
func TestXorCleanOneSq(t *testing.T)   { runXorTables(t, "XorSq", sq, true) }
func TestXorCleanOneNoSq(t *testing.T) { runXorTables(t, "XorNoSq", noSq, true) }

// encryptVals encrypts a full slot vector at the top level.
func encryptVals(t *testing.T, params ckks.Parameters, ecd *ckks.Encoder, enc *rlwe.Encryptor, vals []float64, level int) *rlwe.Ciphertext {
	t.Helper()
	pt := ckks.NewPlaintext(params, level)
	if err := ecd.Encode(vals, pt); err != nil {
		t.Fatalf("encode: %v", err)
	}
	ct, err := enc.EncryptNew(pt)
	if err != nil {
		t.Fatalf("encrypt: %v", err)
	}
	return ct
}
