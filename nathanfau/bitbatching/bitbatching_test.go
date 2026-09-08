package bitbatching

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const testK = 4

// newTestContext builds a standard (complex) CKKS context with relinearization and
// complex-conjugation keys, carrying levels primes above the base one. BitExtract and
// BitExtractInterp need k of them for testK bits, BitExtractClean k+1.
func newTestContext(t *testing.T, levels int) (ckks.Parameters, *ckks.Encoder, *rlwe.Encryptor, *rlwe.Decryptor, *ckks.Evaluator) {
	t.Helper()
	logQ := []int{55}
	for i := 0; i < levels; i++ {
		logQ = append(logQ, 40)
	}
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            logQ,
		LogP:            []int{55, 55},
		LogDefaultScale: 40,
	})
	if err != nil {
		t.Fatalf("params: %v", err)
	}
	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	gk := kgen.GenGaloisKeyNew(params.GaloisElementForComplexConjugation(), sk)
	evk := rlwe.NewMemEvaluationKeySet(rlk, gk)

	ecd := ckks.NewEncoder(params)
	enc := rlwe.NewEncryptor(params, pk)
	dec := rlwe.NewDecryptor(params, sk)
	eval := ckks.NewEvaluator(params, evk)

	debug.EncStd, debug.DecStd, debug.ParamsStd = ecd, dec, params

	return params, ecd, enc, dec, eval
}

// encryptVec encrypts a full complex slot vector.
func encryptVec(t *testing.T, params ckks.Parameters, ecd *ckks.Encoder, enc *rlwe.Encryptor, vals []complex128) *rlwe.Ciphertext {
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

// decodeComplex decrypts ct and returns the complex value of each slot.
func decodeComplex(t *testing.T, ecd *ckks.Encoder, dec *rlwe.Decryptor, ct *rlwe.Ciphertext) []complex128 {
	t.Helper()
	out := make([]complex128, ct.Slots())
	if err := ecd.Decode(dec.DecryptNew(ct), out); err != nil {
		t.Fatalf("decode: %v", err)
	}
	return out
}

// TestBitPack: each of the k planes is a *complex* bit (0/1) + i*(0/1), the real bit from a
// k-bit integer m_A and the imaginary bit from a k-bit integer m_B, as produced by the blockpack
// fold (byte g in the real half, byte g+8 in the imaginary half). Pack the k planes and check
// the packed ciphertext decrypts to m_A + i*m_B, with m_A = sum_i bitA_i*2^i (resp. m_B).
func TestBitPack(t *testing.T) {
	params, ecd, enc, dec, eval := newTestContext(t, 6)
	k := testK
	n := params.MaxSlots()
	rng := rand.New(rand.NewSource(time.Now().UnixNano())) // *random* seed for testing

	mA := make([]int, n)
	mB := make([]int, n)
	planes := make([][]complex128, k)
	for i := range planes {
		planes[i] = make([]complex128, n)
	}
	for s := 0; s < n; s++ {
		mA[s] = rng.Intn(1 << k)
		mB[s] = rng.Intn(1 << k)
		for i := 0; i < k; i++ {
			planes[i][s] = complex(float64((mA[s]>>i)&1), float64((mB[s]>>i)&1))
		}
	}

	for i := 0; i < k; i++ {
		fmt.Printf("planes[%d] = %v\n", i, planes[i][:2])
	}

	bits := make([]*rlwe.Ciphertext, k)
	for i := 0; i < k; i++ {
		bits[i] = encryptVec(t, params, ecd, enc, planes[i])
	}

	for i := 0; i < k; i++ {
		debug.DbgSlotStd("encryptedplanes = ", bits[i])
	}

	ctPack, err := BitPack(eval, bits)
	if err != nil {
		t.Fatalf("bitPack: %v", err)
	}

	debug.DbgSlotStd("BitPack =", ctPack)

	got := decodeComplex(t, ecd, dec, ctPack)
	maxErr := 0.0
	for s := 0; s < n; s++ {
		want := complex(float64(mA[s]), float64(mB[s]))
		if e := cmplx.Abs(got[s] - want); e > maxErr {
			maxErr = e
		}
	}
	t.Logf("bitPack max error (|m_A+i*m_B|) = %.3e", maxErr)
	if maxErr > 0.1 {
		t.Fatalf("bitPack: max error %.3e too large", maxErr)
	}
}

// ---------------------------------------------------------------------------
// extraction
// ---------------------------------------------------------------------------

// bitVariants are the three ways of getting the k bits back out of omega^m, in the order of the
// comparison table: the half-spectrum interpolation already in the pipeline (voie 2), the plain
// interpolation (voie 1 with the cleaning turned off), and the combined interpolation and cleaning
// (voie 1).
var bitVariants = []struct {
	name string
	run  func(ckks.Parameters, *ckks.Evaluator, *rlwe.Ciphertext, int) ([]*rlwe.Ciphertext, error)
}{
	{"voie 2 (BitExtract)", BitExtract},
	{"voie 1 sans clean", BitExtractInterp},
	{"voie 1 avec clean", BitExtractClean},
}

// extractExps are the input error sizes the table sweeps, one row each. They stop at 14 because
// past that the quadratic output would sit under the noise floor of a fresh ciphertext at scale
// 2^40 and the two regimes could no longer be told apart, and they start at 6 because h is only
// defined for |eps| <= 1/(t-1), i.e. 2^-3.9 at k = 4.
//
// Every row carries an injected error, on purpose: on an EXACT input the three variants are
// indistinguishable, there being nothing left for the cleaning to clean.
var extractExps = []int{6, 10, 14}

// encodeRoots returns omega^m for a random m per slot, perturbed by an error of modulus 2^-exp in a
// random direction (exp <= 0 leaves them exact), the EXACT roots to measure that input against, and
// the bit planes they encode.
func encodeRoots(rng *rand.Rand, n, k, exp int) (vals, exact []complex128, planes [][]float64) {
	t := 1 << k
	base := 2 * math.Pi / float64(t)
	vals = make([]complex128, n)
	exact = make([]complex128, n)
	planes = make([][]float64, k)
	for i := range planes {
		planes[i] = make([]float64, n)
	}
	for s := 0; s < n; s++ {
		m := rng.Intn(t)
		exact[s] = cmplx.Exp(complex(0, base*float64(m)))
		vals[s] = exact[s]
		if exp > 0 {
			vals[s] += cmplx.Rect(math.Pow(2, -float64(exp)), rng.Float64()*2*math.Pi)
		}
		for i := 0; i < k; i++ {
			planes[i][s] = float64((m >> i) & 1)
		}
	}
	return vals, exact, planes
}

// TestBitExtract runs the three extractions on the SAME ciphertext, at several input error sizes,
// and prints precision, levels and time side by side.
//
// Every variant has to return every bit right, on every row -- that is the extraction working at
// all. The table is what says whether voie 1 is worth its extra prime:
//
//   - voie 2 and voie 1 sans clean are both LINEAR, so their precision tracks the input one for one;
//   - voie 1 avec clean is QUADRATIC, so its lead over the other two doubles every time the input
//     gets two bits cleaner, until the ciphertext noise floor stops it.
func TestBitExtract(t *testing.T) {
	params, ecd, enc, dec, eval := newTestContext(t, 8)
	k, n := testK, params.MaxSlots()
	rng := rand.New(rand.NewSource(time.Now().UnixNano())) // *random* seed for testing

	groups := make([]utils.TableGroup, len(bitVariants))
	for i, v := range bitVariants {
		groups[i] = utils.TableGroup{Name: v.name, Cols: []string{"avg", "worst", "lv", "time"}}
	}
	leads := make([]string, len(extractExps))
	rows := make([][]string, len(extractExps))
	ahead := make([]float64, len(extractExps)) // voie 1 avec clean, minus voie 2

	for ei, exp := range extractExps {
		vals, exact, planes := encodeRoots(rng, n, k, exp)
		ct := encryptVec(t, params, ecd, enc, vals)

		// Against the EXACT roots, not the jittered plaintext: the injected error is what the
		// variants amplify, so it is what the input column has to show.
		in, err := utils.BitDistanceCt(ecd, dec, ct, exact, 0)
		if err != nil {
			t.Fatalf("exp=%d input precision: %v", exp, err)
		}
		leads[ei] = fmt.Sprintf("2^-%d (%.1f/%.1f)", exp, in.AvgPrec, -math.Log2(in.Worst))

		fmt.Printf("vals = %f\n", vals[:2])
		debug.DbgSlotStd("omega^m =", ct)
		debug.DbgChain("Chain before extraction :", eval, ct)

		row := make([]string, 0, 4*len(bitVariants))
		avg := make([]float64, len(bitVariants))
		for vi, v := range bitVariants {
			t0 := time.Now()
			bits, err := v.run(params, eval, ct.CopyNew(), k)
			if err != nil {
				t.Fatalf("exp=%d %s: %v", exp, v.name, err)
			}
			dur := time.Since(t0)
			if len(bits) != k {
				t.Fatalf("exp=%d %s: returned %d bits, want %d", exp, v.name, len(bits), k)
			}

			// The labels carry the variant, the three being interleaved in one run.
			for i := 0; i < k; i++ {
				debug.DbgSlotStd("bitsExtracted "+v.name+" = ", bits[i])
			}
			debug.DbgChain("Chain after extraction "+v.name+" :", eval, bits[0])

			sumAvg, worst := 0.0, 0.0
			for i := 0; i < k; i++ {
				// All three are contracted to come out on the input's scale. A drift there would be
				// a silent precision floor and nothing else, since ckks.Add takes the INTEGER part
				// of a scale ratio and so does nothing for a ratio of 1.0002.
				if d := math.Abs(bits[i].Scale.Float64()/ct.Scale.Float64() - 1); d > math.Pow(2, -40) {
					t.Errorf("exp=%d %s bit %d: output scale is 2^%.2f off the input scale, the bookkeeping drifted",
						exp, v.name, i, math.Log2(d))
				}
				st, err := utils.BitDistanceCt(ecd, dec, bits[i], planes[i], 0.5)
				if err != nil {
					t.Fatalf("exp=%d %s bit %d: %v", exp, v.name, i, err)
				}
				if st.Wrong != 0 {
					t.Errorf("exp=%d %s bit %d: %d slots farther than 0.5 from their value (worst %.4f at slot %d)",
						exp, v.name, i, st.Wrong, st.Worst, st.WorstSlot)
				}
				sumAvg += st.AvgPrec
				worst = math.Max(worst, st.Worst)
			}
			avg[vi] = sumAvg / float64(k)
			row = append(row,
				fmt.Sprintf("%.2f", avg[vi]), fmt.Sprintf("%.2f", -math.Log2(worst)),
				fmt.Sprintf("%d->%d", ct.Level(), bits[0].Level()),
				dur.Round(time.Millisecond).String())
		}
		rows[ei] = row
		ahead[ei] = avg[2] - avg[0]
	}

	fmt.Printf("Extraction de %d bits, logN %d, %d slots\n%s",
		k, params.LogN(), n, utils.BoxTable("entree (avg/worst)", groups, leads, rows))
	fmt.Printf("avance de « %s » sur « %s » : %s\n", bitVariants[2].name, bitVariants[0].name, fmtAhead(ahead))

	// The cleaning has to buy something at every size of the sweep.
	for ei, exp := range extractExps {
		if ahead[ei] < 2 {
			t.Errorf("exp=%d : « %s » ne devance « %s » que de %.2f bits, attendu >= 2",
				exp, bitVariants[2].name, bitVariants[0].name, ahead[ei])
		}
	}
	// And the lead has to GROW as the input gets cleaner: one bit of lead per bit of input, which is
	// the quadratic law and the whole justification for voie 1's extra prime. Anything under half of
	// the span means the cleaning is not quadratic -- the article's t/2 threshold, typically.
	last := len(extractExps) - 1
	span := float64(extractExps[last] - extractExps[0])
	if gain := ahead[last] - ahead[0]; gain < span/2 {
		t.Errorf("l'avance n'a gagne que %.2f bits sur %.0f bits d'entree, attendu >= %.0f : l'erreur n'est pas quadratique",
			gain, span, span/2)
	}
}

func fmtAhead(v []float64) string {
	s := ""
	for i, x := range v {
		if i > 0 {
			s += ", "
		}
		s += fmt.Sprintf("2^-%d: %+.2f", extractExps[i], x)
	}
	return s
}
