package aes_test

// The SubBytes variants against their exactly aligned counterparts (sbox_exact_he.go), in the
// conditions of the transciphering pipeline: its chain (params2), the conjugate-invariant ring, the
// level SubBytes starts at. Only the relinearization key is generated, so logN 16 fits on a machine
// that could not bootstrap.
//
// TestSubBytesExact runs them on a fresh state; TestSubBytesExactNoise on the same state carrying
// an input error of 2^-e, e swept by -noise, the way the cleaning test and nathanfau/bench inject it.
//
//	go test ./nathanfau/aes/ -run '^TestSubBytesExactNoise$' -v -timeout 0
//	go test ./nathanfau/aes/ -run '^TestSubBytesExactNoise$' -v -timeout 0 -noise 0,20,12,8 -sbv 2,3
//	go test ./nathanfau/aes/ -run '^TestSubBytesExact$' -v -timeout 0 -logn 16 -sbv 3
//
// The floor the exact version removes grows with logN: the primes, = 1 mod 4N, drift further from
// their power of two.

import (
	"flag"
	"fmt"
	"math"
	"math/rand"
	"strconv"
	"strings"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockpack"
	"github.com/tuneinsight/lattigo/v6/nathanfau/params2"
	"github.com/tuneinsight/lattigo/v6/nathanfau/transciphering"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var (
	exLogN  = flag.Int("logn", 12, "ring degree")
	exLogQi = flag.Int("logqi", params2.LogScale, "size of the chain primes and of the scale")
	exSTC   = flag.String("stc", "2x30", `SlotsToCoeffs primes, "2x30" or "30,30"; it moves the level SubBytes starts at`)
	exSbv   = flag.String("sbv", "2", "SubBytes versions to compare, among 1,2,3 (1 is the slowest)")
	exSeed  = flag.Int64("seed", 42, "seed of the block draw and of the injected errors")
	exNoise = flag.String("noise", "0,30,26,22,18,14,12,10,8,6", "TestSubBytesExactNoise: input error sizes 2^-e; 0 is a fresh ciphertext")
)

// exactEnv is everything both tests share: the pipeline's chain in the CI ring, the keys, the
// batch and its expected SubBytes.
type exactEnv struct {
	params    ckks.Parameters
	logSTC    []int
	level     int
	ecd       *ckks.Encoder
	enc       *rlwe.Encryptor
	dec       *rlwe.Decryptor
	ae        *aes.Evaluator
	versions  []int
	blocks    [][16]byte
	want      [][16]byte  // SubBytes of every block
	inSlots   [][]float64 // the 64 input bit planes
	wantSlots [][]float64 // the 64 expected output bit planes
}

func newExactEnv(t *testing.T) *exactEnv {
	t.Helper()
	if testing.Short() {
		t.Skip("SubBytes on a full batch at the pipeline's ring degree")
	}
	e := &exactEnv{}
	var err error
	if e.logSTC, err = params2.ParseLogSTC(*exSTC); err != nil {
		t.Fatalf("-stc: %v", err)
	}
	for _, f := range strings.Split(*exSbv, ",") {
		v, err := strconv.Atoi(strings.TrimSpace(f))
		if err != nil || v < 1 || v > 3 {
			t.Fatalf("-sbv: %q, want versions among 1,2,3", f)
		}
		e.versions = append(e.versions, v)
	}

	// The pipeline's chain, in the CI ring SubBytes runs in (built the way convctx builds it).
	std, _, err := params2.TranscipheringParamsWith(*exLogN, 4, 2, 5, params2.Shape{LogSTC: e.logSTC, LogQi: *exLogQi})
	if err != nil {
		t.Fatalf("params: %v", err)
	}
	lit := std.ParametersLiteral()
	lit.RingType = ring.ConjugateInvariant
	if e.params, err = ckks.NewParametersFromLiteral(lit); err != nil {
		t.Fatalf("CI params: %v", err)
	}
	// SubBytes starts where the pipeline hands it the state: SubBytesLevel, lifted by every
	// SlotsToCoeffs prime beyond the first (transciphering.NewContextWith2 does the same).
	e.level = transciphering.SubBytesLevel + len(e.logSTC) - 1

	kgen := rlwe.NewKeyGenerator(e.params)
	sk := kgen.GenSecretKeyNew()
	eval := ckks.NewEvaluator(e.params, rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk)))
	e.ecd = ckks.NewEncoder(e.params)
	e.enc = rlwe.NewEncryptor(e.params, sk)
	e.dec = rlwe.NewDecryptor(e.params, sk)
	e.ae = aes.NewEvaluator(eval)

	e.blocks = blockpack.RandomBlocks(e.params, rand.New(rand.NewSource(*exSeed)))
	e.want = make([][16]byte, len(e.blocks))
	for i, b := range e.blocks {
		e.want[i] = b
		aes.SubBytes(e.want[i][:])
	}
	e.inSlots = blockpack.WantSlots(e.params, e.blocks)
	e.wantSlots = blockpack.WantSlots(e.params, e.want)

	fmt.Printf("chain: logN=%d, q_i=%d, STC=%s, SubBytes input level %d, %d blocks, seed %d\n",
		*exLogN, *exLogQi, params2.FormatLogSTC(e.logSTC), e.level, len(e.blocks), *exSeed)
	fmt.Printf("  log2 q_i: %v\n", e.params.LogQi())
	return e
}

// encrypt packs the batch at SubBytes' input level, every slot of every bit plane moved by an error
// of the order of 2^-exp (0: none), with the same law as the cleaning test: a random sign and a
// magnitude jittered by up to min(exp/32, 0.9).
func (e *exactEnv) encrypt(t *testing.T, exp int, rng *rand.Rand) blockpack.Packed {
	t.Helper()
	mag, jitter := math.Exp2(-float64(exp)), math.Min(float64(exp)/32, 0.9)
	var p blockpack.Packed
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			vals := append([]float64(nil), e.inSlots[8*g+b]...)
			if exp > 0 {
				for s := range vals {
					sign := float64(rng.Intn(2)*2 - 1)
					vals[s] += sign * mag * (1 + (rng.Float64()*2-1)*jitter)
				}
			}
			pt := ckks.NewPlaintext(e.params, e.level)
			if err := e.ecd.Encode(vals, pt); err != nil {
				t.Fatalf("encode: %v", err)
			}
			ct := ckks.NewCiphertext(e.params, 1, e.level)
			if err := e.enc.Encrypt(pt, ct); err != nil {
				t.Fatalf("encrypt: %v", err)
			}
			p[g][b] = ct
		}
	}
	return p
}

// sboxRun is one circuit on the 8 groups of one state.
type sboxRun struct {
	dur         time.Duration
	relin, resc int // per byte, the circuit's own
	align       int // per byte, the rescales spent aligning scales on level drops
	levels      string
	scales      int     // distinct scales among the 64 outputs
	worst, avg  float64 // bits
	wrong       int     // blocks the cleartext S-box disagrees with
}

func (e *exactEnv) run(t *testing.T, name string, fn func(aes.ByteHE) (aes.ByteHE, error), st blockpack.Packed) sboxRun {
	t.Helper()
	utils.ResetOps()
	var out blockpack.Packed
	var err error
	t0 := time.Now()
	for g := range 8 {
		if out[g], err = fn(st[g]); err != nil {
			t.Fatalf("%s group %d: %v", name, g, err)
		}
	}
	r := sboxRun{dur: time.Since(t0), relin: utils.Ops.Relin / 8, resc: utils.Ops.Rescale / 8, align: utils.Ops.Align / 8}

	lvMin, lvMax := out[0][0].Level(), out[0][0].Level()
	scales := map[string]bool{}
	for _, ct := range out.Cts() {
		lvMin, lvMax = min(lvMin, ct.Level()), max(lvMax, ct.Level())
		scales[ct.Scale.Value.Text('g', 40)] = true
	}
	r.levels, r.scales = strconv.Itoa(lvMin), len(scales)
	if lvMax != lvMin {
		r.levels = fmt.Sprintf("%d-%d", lvMin, lvMax)
	}

	stats, err := utils.BitDistance(e.ecd, e.dec, out.Cts(), e.wantSlots, 0)
	if err != nil {
		t.Fatalf("%s: %v", name, err)
	}
	agg, _ := utils.WorstOf(stats)
	r.worst, r.avg = -math.Log2(agg.Worst), agg.AvgPrec

	got, _, err := blockpack.Decrypt(e.params, e.ecd, e.dec, out)
	if err != nil {
		t.Fatalf("%s decrypt: %v", name, err)
	}
	for i := range e.want {
		if got[i] != e.want[i] {
			r.wrong++
		}
	}
	return r
}

// inputPrec is the precision of a packed input against its exact bits, in bits (worst, mean).
func (e *exactEnv) inputPrec(t *testing.T, st blockpack.Packed) (worst, avg float64) {
	t.Helper()
	stats, err := utils.BitDistance(e.ecd, e.dec, st.Cts(), e.inSlots, 0)
	if err != nil {
		t.Fatalf("input precision: %v", err)
	}
	agg, _ := utils.WorstOf(stats)
	return -math.Log2(agg.Worst), agg.AvgPrec
}

func (e *exactEnv) original(v int) func(aes.ByteHE) (aes.ByteHE, error) {
	return func(b aes.ByteHE) (aes.ByteHE, error) { return e.ae.SubByte(b, v) }
}

func (e *exactEnv) exact(v int) func(aes.ByteHE) (aes.ByteHE, error) {
	return func(b aes.ByteHE) (aes.ByteHE, error) { return e.ae.SubByteExact(b, v) }
}

// TestSubBytesExact compares each SubBytes variant with its exactly aligned counterpart on the SAME
// fresh state.
func TestSubBytesExact(t *testing.T) {
	e := newExactEnv(t)
	st := e.encrypt(t, 0, nil)

	var leads []string
	var rows [][]string
	for _, v := range e.versions {
		for _, c := range []struct {
			name string
			fn   func(aes.ByteHE) (aes.ByteHE, error)
		}{{fmt.Sprintf("V%d", v), e.original(v)}, {fmt.Sprintf("V%d exact", v), e.exact(v)}} {
			r := e.run(t, c.name, c.fn, st)
			if r.wrong != 0 {
				t.Errorf("%s: %d blocks wrong on a fresh state", c.name, r.wrong)
			}
			fmt.Printf("  %-9s done in %s, worst %.2f bits\n", c.name, r.dur.Round(time.Millisecond), r.worst)
			leads = append(leads, c.name)
			rows = append(rows, []string{
				strconv.Itoa(r.relin), strconv.Itoa(r.resc), strconv.Itoa(r.align), r.dur.Round(time.Millisecond).String(),
				r.levels, strconv.Itoa(r.scales), f2(r.worst), f2(r.avg), strconv.Itoa(r.wrong),
			})
		}
	}
	groups := []utils.TableGroup{
		{Name: "ops / octet", Cols: []string{"relin", "rescale", "align"}},
		{Name: "8 octets", Cols: []string{"temps"}},
		{Name: "sortie", Cols: []string{"niveau", "#echelles"}},
		{Name: "precision (bits)", Cols: []string{"pire", "moy"}},
		{Name: "blocs", Cols: []string{"faux"}},
	}
	fmt.Printf("SubBytes, alignement exact des echelles\n%s", utils.BoxTable("circuit", groups, leads, rows))
	fmt.Println(`  rescale / octet : les rescales du circuit ; align : les rescales d'alignement de la version exacte.
  #echelles : nombre d'echelles distinctes parmi les 64 sorties ; 1 = toutes identiques.`)
}

// TestSubBytesExactNoise runs each version and its exact counterpart on the SAME state for every
// input error size of -noise, and reports time, precision and wrong blocks side by side. Nothing
// fails on wrong blocks: past some error that is what is being measured.
func TestSubBytesExactNoise(t *testing.T) {
	e := newExactEnv(t)
	var exps []int
	for _, f := range strings.Split(*exNoise, ",") {
		x, err := strconv.Atoi(strings.TrimSpace(f))
		if err != nil || x < 0 {
			t.Fatalf("-noise: %q, want exponents e >= 0 (error 2^-e, 0 = fresh)", f)
		}
		exps = append(exps, x)
	}

	for _, v := range e.versions {
		orig, exact := fmt.Sprintf("V%d", v), fmt.Sprintf("V%d exact", v)
		var leads []string
		var rows [][]string
		for _, x := range exps {
			label := "frais"
			if x > 0 {
				label = fmt.Sprintf("2^-%d", x)
			}
			st := e.encrypt(t, x, rand.New(rand.NewSource(*exSeed+int64(x))))
			inWorst, _ := e.inputPrec(t, st)
			a := e.run(t, orig, e.original(v), st)
			b := e.run(t, exact, e.exact(v), st)
			fmt.Printf("  %-6s entree %5.2f | %s %s pire %5.2f | %s %s pire %5.2f\n", label, inWorst,
				orig, a.dur.Round(time.Millisecond), a.worst, exact, b.dur.Round(time.Millisecond), b.worst)

			leads = append(leads, label)
			rows = append(rows, []string{
				f2(inWorst),
				a.dur.Round(time.Millisecond).String(), f2(a.worst), f2(a.avg), strconv.Itoa(a.wrong),
				b.dur.Round(time.Millisecond).String(), f2(b.worst), f2(b.avg), strconv.Itoa(b.wrong),
				fmt.Sprintf("%+.2f", b.worst-a.worst), fmt.Sprintf("%+.2f", b.avg-a.avg),
				fmt.Sprintf("%+.0f%%", 100*(float64(b.dur)/float64(a.dur)-1)),
			})
		}
		groups := []utils.TableGroup{
			{Name: "entree", Cols: []string{"pire"}},
			{Name: orig, Cols: []string{"temps", "pire", "moy", "faux"}},
			{Name: exact, Cols: []string{"temps", "pire", "moy", "faux"}},
			{Name: "exact - original", Cols: []string{"pire", "moy", "temps"}},
		}
		fmt.Printf("\nSubBytes %s contre %s, meme chiffre, erreur d'entree 2^-e\n%s", orig, exact,
			utils.BoxTable("erreur", groups, leads, rows))
	}
	fmt.Println(`  precisions en bits, contre le SubBytes en clair ; entree = precision de l'etat avant SubBytes.
  faux = blocs (sur toute la batch) dont au moins un bit sort du mauvais cote de 1/2.`)
}

func f2(x float64) string {
	if math.IsInf(x, 0) || math.IsNaN(x) {
		return "-"
	}
	return fmt.Sprintf("%.2f", x)
}
