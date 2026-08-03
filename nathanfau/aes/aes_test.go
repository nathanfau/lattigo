package aes

//	go test ./nathanfau/aes/ -v
//
// We deliberately do not test ShiftRows: it consumes the nibbles ALREADY split into their Re / Im
// streams by the Algo1 pause, so exercising it here would mean reproducing that split outside the
// refresh it belongs to.

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockpack"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// testBlocks fills the whole batch with independent random AES blocks, one per slot,
//
//	and logs the seed so a failure can be replayed.
func testBlocks(t *testing.T, params ckks.Parameters) [][16]byte {
	t.Helper()
	seed := time.Now().UnixNano()
	t.Logf("random seed = %d", seed)
	return blockpack.RandomBlocks(params, rand.New(rand.NewSource(seed)))
}

func newAESContext(t *testing.T) (ckks.Parameters, *ckks.Encoder, *rlwe.Encryptor, *rlwe.Decryptor, *Evaluator) {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            12,
		LogQ:            []int{55, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{61},
		LogDefaultScale: 45,
	})
	if err != nil {
		t.Fatalf("params: %v", err)
	}
	debug.DbgParams("aes", params)
	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	evk := rlwe.NewMemEvaluationKeySet(rlk)

	ecd := ckks.NewEncoder(params)
	enc := rlwe.NewEncryptor(params, pk)
	dec := rlwe.NewDecryptor(params, sk)
	ae := NewEvaluator(ckks.NewEvaluator(params, evk))

	debug.EncStd, debug.DecStd, debug.ParamsStd = ecd, dec, params

	return params, ecd, enc, dec, ae
}

// worstPrec is the precision of the worst slot of the whole packed state, in bits.
func worstPrec(t *testing.T, ecd *ckks.Encoder, dec *rlwe.Decryptor, p blockpack.Packed, want [][]float64) float64 {
	t.Helper()
	st, err := utils.BitDistance(ecd, dec, p.Cts(), want, 0)
	if err != nil {
		t.Fatalf("BitDistance: %v", err)
	}
	agg, _ := utils.WorstOf(st)
	return -math.Log2(agg.Worst)
}

// subByteVariant is one SubBytes implementation under test, with its expected
// relin==rescale count per S-box circuit (the BCKK op-count baseline).
type subByteVariant struct {
	name        string
	fn          func(ByteHE) (ByteHE, error)
	wantPerByte int
}

// TestSubBytes runs every SubBytes variant on the same batch of blocks as sub-tests
func TestSubBytes(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)

	blocks := testBlocks(t, params)
	st, err := blockpack.Encrypt(params, ecd, enc, blocks, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc state: %v", err)
	}
	want := make([][16]byte, len(blocks))
	for i, blk := range blocks {
		want[i] = blk
		SubBytes(want[i][:])
	}

	variants := []subByteVariant{
		{"V1", ae.SubByteV1, 247}, // non-lazy: one relin+rescale per degree-2+ monomial (255-8)
		{"V2", ae.SubByteV2, 98},  // lazy deg>=4: 61 factors + 29 low-degree leaves + 8 accumulators
		{"V3", ae.SubByteV3, 69},  // lazy: 61 reused factors + 8 output-bit accumulators
	}

	var leads []string
	var rows [][]string
	for _, v := range variants {
		t.Run(v.name, func(t *testing.T) {
			utils.ResetOps()
			var out blockpack.Packed
			t0 := time.Now()
			for g := range 8 {
				var e error
				if out[g], e = v.fn(st[g]); e != nil {
					t.Fatalf("%s group %d: %v", v.name, g, e)
				}
			}
			dur := time.Since(t0)
			ops := utils.Ops

			if e := blockpack.Check(params, ecd, dec, out, want); e != nil {
				t.Fatalf("%s %v", v.name, e)
			}
			if ops.Relin != 8*v.wantPerByte || ops.Rescale != 8*v.wantPerByte {
				t.Errorf("%s op count = %d relin / %d rescale, want %d / %d (%d per circuit)",
					v.name, ops.Relin, ops.Rescale, 8*v.wantPerByte, 8*v.wantPerByte, v.wantPerByte)
			}

			leads = append(leads, v.name)
			rows = append(rows, []string{
				fmt.Sprintf("%d", ops.Relin/8), fmt.Sprintf("%d", ops.Relin),
				fmt.Sprintf("%d", ops.Rescale/8), fmt.Sprintf("%d", ops.Rescale),
				dur.Round(time.Millisecond).String(),
				fmt.Sprintf("%.2f", worstPrec(t, ecd, dec, out, blockpack.WantSlots(params, want))),
			})
		})
	}

	groups := []utils.TableGroup{
		{Name: "relin", Cols: []string{"/circ", "total"}},
		{Name: "rescale", Cols: []string{"/circ", "total"}},
		{Name: "circuit", Cols: []string{"time", "worst prec(b)"}},
	}
	fmt.Printf("SubBytes variants\n%s", utils.BoxTable("SubBytes", groups, leads, rows))
}

func TestMixColumns(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)

	blocks := testBlocks(t, params)
	st, err := blockpack.Encrypt(params, ecd, enc, blocks, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc state: %v", err)
	}

	want := make([][16]byte, len(blocks))
	for i, blk := range blocks {
		want[i] = MixColumnsRM(blk)
	}

	debug.DbgBitStd("in bit[0][0]", ae.eval, st[0][0], blockpack.SlotVec(params, blocks, 0, 0))

	out, err := ae.MixColumns(st)
	if err != nil {
		t.Fatalf("MixColumns: %v", err)
	}

	debug.DbgBitStd("out bit[0][0]", ae.eval, out[0][0], blockpack.SlotVec(params, want, 0, 0))

	if err := blockpack.Check(params, ecd, dec, out, want); err != nil {
		t.Fatalf("MixColumns %v", err)
	}
}

// xorVariant is one XOR circuit under test; pick binds it to an evaluator.
type xorVariant struct {
	name string
	pick func(*Evaluator) XorFunc
}

var xorVariants = []xorVariant{
	{"(x-y)^2", func(a *Evaluator) XorFunc { return a.XorSq }},
	{"x+y-2xy", func(a *Evaluator) XorFunc { return a.XorNoSq }},
}

// Input error sizes: the same ladder as the cleaning package's sweep, so the tables read against
// each other. The S-box ladder starts later -- its ANF is degree 7, so it blows a coarse input up
// well past the flip threshold and the row would say nothing.
var (
	xorExps = []int{2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30}
	mixExps = []int{26, 20, 16, 14, 12}
)

// jitterNoise adds an error of magnitude about 2^-exp to every value, sign drawn at random.
func jitterNoise(rng *rand.Rand, vals []float64, exp int) {
	mag := math.Exp2(-float64(exp))
	frac := math.Min(0.5, float64(exp)/32)
	for i := range vals {
		sign := float64(rng.Intn(2)*2 - 1)
		vals[i] += sign * mag * (1 + (rng.Float64()*2-1)*frac)
	}
}

// encVals encrypts a full slot vector at the top level.
func encVals(t *testing.T, params ckks.Parameters, ecd *ckks.Encoder, enc *rlwe.Encryptor, vals []float64) *rlwe.Ciphertext {
	t.Helper()
	pt := ckks.NewPlaintext(params, params.MaxLevel())
	if err := ecd.Encode(vals, pt); err != nil {
		t.Fatalf("encode: %v", err)
	}
	ct, err := enc.EncryptNew(pt)
	if err != nil {
		t.Fatalf("encrypt: %v", err)
	}
	return ct
}

// TestXorVariants checks the truth table and the level cost of both XOR circuits, then sweeps the
// input error and compares their precision and time side by side.
func TestXorVariants(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)
	n := params.MaxSlots()
	seed := time.Now().UnixNano()
	t.Logf("random seed = %d", seed)

	// Truth table and level cost, on inputs carrying nothing but the encryption noise.
	rng := rand.New(rand.NewSource(seed))
	bx, by := make([]float64, n), make([]float64, n)
	want := make([]float64, n)
	for s := range want {
		bx[s], by[s] = float64(rng.Intn(2)), float64(rng.Intn(2))
		want[s] = float64(int(bx[s]) ^ int(by[s]))
	}
	cx := encVals(t, params, ecd, enc, bx)
	cy := encVals(t, params, ecd, enc, by)

	for _, v := range xorVariants {
		out, err := v.pick(ae)(cx, cy)
		if err != nil {
			t.Fatalf("%s: %v", v.name, err)
		}
		if got := params.MaxLevel() - out.Level(); got != 1 {
			t.Errorf("%s consumed %d levels, want 1", v.name, got)
		}
		st, err := utils.BitDistanceCt(ecd, dec, out, want, 0.25)
		if err != nil {
			t.Fatalf("%s: %v", v.name, err)
		}
		if st.Wrong != 0 {
			t.Errorf("%s: %d/%d slots off the XOR truth table by more than 0.25", v.name, st.Wrong, st.Slots)
		}
	}

	// Error sweep. Both variants see the SAME noisy pair at each size, so a row compares straight
	// across and the input reference is measured once.
	leads := make([]string, len(xorExps))
	rows := make([][]string, len(xorExps))

	for ei, e := range xorExps {
		rng := rand.New(rand.NewSource(int64(e) + 1))
		bx, by := make([]float64, n), make([]float64, n)
		nx, ny := make([]float64, n), make([]float64, n)
		want := make([]float64, n)
		for s := range want {
			bx[s], by[s] = float64(rng.Intn(2)), float64(rng.Intn(2))
			want[s] = float64(int(bx[s]) ^ int(by[s]))
		}
		copy(nx, bx)
		copy(ny, by)
		jitterNoise(rng, nx, e)
		jitterNoise(rng, ny, e)
		cx := encVals(t, params, ecd, enc, nx)
		cy := encVals(t, params, ecd, enc, ny)

		// Reference: the worse of the two input bits.
		sx, err := utils.BitDistanceCt(ecd, dec, cx, bx, 0)
		if err != nil {
			t.Fatalf("input x (%s): %v", fmt.Sprintf("2^-%d", e), err)
		}
		sy, err := utils.BitDistanceCt(ecd, dec, cy, by, 0)
		if err != nil {
			t.Fatalf("input y (%s): %v", fmt.Sprintf("2^-%d", e), err)
		}
		inAvg := math.Min(sx.AvgPrec, sy.AvgPrec)
		inWorst := -math.Log2(math.Max(sx.Worst, sy.Worst))
		leads[ei] = fmt.Sprintf("2^-%d %.1f/%.1f", e, inAvg, inWorst)

		row := make([]string, 0, 5*len(xorVariants))
		for _, v := range xorVariants {
			xor := v.pick(ae)
			t0 := time.Now()
			out, err := xor(cx, cy)
			if err != nil {
				t.Fatalf("%s (2^-%d): %v", v.name, e, err)
			}
			dur := time.Since(t0)
			st, err := utils.BitDistanceCt(ecd, dec, out, want, 0)
			if err != nil {
				t.Fatalf("%s (%s): %v", v.name, fmt.Sprintf("2^-%d", e), err)
			}
			avg, worst := st.AvgPrec, -math.Log2(st.Worst)
			row = append(row,
				fmt.Sprintf("%.2f", avg), fmt.Sprintf("%+.2f", avg-inAvg),
				fmt.Sprintf("%.2f", worst), fmt.Sprintf("%+.2f", worst-inWorst),
				dur.Round(time.Microsecond).String())
		}
		rows[ei] = row
	}

	groups := make([]utils.TableGroup, len(xorVariants))
	for i, v := range xorVariants {
		groups[i] = utils.TableGroup{
			Name: v.name,
			Cols: []string{"avg", "diff avg", "worst", "diff worst", "time"},
		}
	}
	fmt.Printf("XOR variants, one gate\n%s",
		utils.BoxTable("in / avg/worst", groups, leads, rows))
}

// noisyPack encrypts the bit-sliced state of blocks with an error of about 2^-exp on every slot.
func noisyPack(t *testing.T, params ckks.Parameters, ecd *ckks.Encoder, enc *rlwe.Encryptor,
	blocks [][16]byte, exp int, rng *rand.Rand) blockpack.Packed {
	t.Helper()
	var p blockpack.Packed
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			vals := blockpack.SlotVec(params, blocks, g, b)
			jitterNoise(rng, vals, exp)
			p[g][b] = encVals(t, params, ecd, enc, vals)
		}
	}
	return p
}

// TestMixColumnsXor compares the two XOR circuits where the pipeline actually spends them: the
// balanced trees of MixColumns, three levels deep over 5 to 7 leaves. Both are checked against
// FIPS-197 on exact inputs, then measured as the input error grows.
//
// "bad" counts the slots that would round to the wrong bit, over the whole 64-ciphertext state.
func TestMixColumnsXor(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)
	blocks := testBlocks(t, params)

	want := make([][16]byte, len(blocks))
	for i, blk := range blocks {
		want[i] = MixColumnsRM(blk)
	}
	wantSlots := blockpack.WantSlots(params, want)

	leads := make([]string, len(mixExps))
	rows := make([][]string, len(mixExps))

	for ei, e := range mixExps {
		st := noisyPack(t, params, ecd, enc, blocks, e, rand.New(rand.NewSource(int64(e)+1)))
		leads[ei] = fmt.Sprintf("2^-%d", e)

		row := make([]string, 0, 3*len(xorVariants))
		for _, v := range xorVariants {
			t0 := time.Now()
			out, err := ae.MixColumnsWith(v.pick(ae), st)
			dur := time.Since(t0)
			if err != nil {
				t.Fatalf("%s (%s): %v", v.name, fmt.Sprintf("2^-%d", e), err)
			}
			if got := params.MaxLevel() - out[0][0].Level(); got != 3 {
				t.Errorf("%s: MixColumns consumed %d levels, want 3", v.name, got)
			}

			// tol 0.25: slots that would round to the wrong bit are a failure, not a column.
			all, err := utils.BitDistance(ecd, dec, blockpack.Packed(out).Cts(), wantSlots, 0.25)
			if err != nil {
				t.Fatalf("%s (%s): %v", v.name, fmt.Sprintf("2^-%d", e), err)
			}
			agg, _ := utils.WorstOf(all)
			for _, s := range all {
				if s.Wrong != 0 {
					t.Errorf("%s (%s): %d/%d slots off the MixColumns oracle", v.name, fmt.Sprintf("2^-%d", e), s.Wrong, s.Slots)
					break
				}
			}
			row = append(row,
				fmt.Sprintf("%.2f", agg.AvgPrec), fmt.Sprintf("%.2f", -math.Log2(agg.Worst)),
				dur.Round(time.Millisecond).String())
		}
		rows[ei] = row
	}

	groups := make([]utils.TableGroup, len(xorVariants))
	for i, v := range xorVariants {
		groups[i] = utils.TableGroup{
			Name: v.name,
			Cols: []string{"avg", "worst", "time"},
		}
	}
	fmt.Printf("MixColumns, XOR trees over %d slots x 64 ct\n%s",
		params.MaxSlots(), utils.BoxTable("in err", groups, leads, rows))
}
