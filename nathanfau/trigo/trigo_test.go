package trigo

import (
	"math"
	"math/cmplx"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// test parameters shared by every trigo test.
const (
	testK      = 1   // Chebyshev interval [-K, K]
	testDegree = 31  // Chebyshev degree
	testR      = 3   // number of double-angle squarings
	testPeriod = 1.0 // target frequency is 1/period
)

func newTrigoContext(t *testing.T) (ckks.Parameters, *ckks.Encoder, *rlwe.Encryptor, *rlwe.Decryptor, *ckks.Evaluator) {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            14,
		LogQ:            []int{50, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{61},
		LogDefaultScale: 45,
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

// encryptInput encrypts n real values drawn uniformly in [-K, K] (fixed seed) and
// returns the ciphertext together with the plaintext values.
func encryptInput(t *testing.T, params ckks.Parameters, ecd *ckks.Encoder, enc *rlwe.Encryptor) (*rlwe.Ciphertext, []float64) {
	t.Helper()
	n := params.MaxSlots()
	rng := rand.New(rand.NewSource(time.Now().UnixNano())) // *random* seed for testing
	x := make([]float64, n)
	for s := 0; s < n; s++ {
		x[s] = (rng.Float64()*2 - 1) * float64(testK)
	}
	pt := ckks.NewPlaintext(params, params.MaxLevel())
	if err := ecd.Encode(x, pt); err != nil {
		t.Fatalf("encode: %v", err)
	}
	ct, err := enc.EncryptNew(pt)
	if err != nil {
		t.Fatalf("encrypt: %v", err)
	}
	return ct, x
}

func TestEvalCos(t *testing.T) {
	params, ecd, enc, dec, eval := newTrigoContext(t)
	ct, x := encryptInput(t, params, ecd, enc)

	fT := testPeriod * math.Pow(2, float64(testR))
	want := make([]float64, len(x))
	for s := range x {
		want[s] = math.Cos(2 * math.Pi * x[s] / fT)
	}

	debug.DbgSlotStd("input x =", ct)
	debug.DbgChain("Chain before EvalCos :", eval, ct)
	debug.PrintPrecStd("Prec in: ", x, ct)
	out, err := EvalCos(params, eval, ct, testK, testPeriod, testR, testDegree)
	if err != nil {
		t.Fatalf("EvalCos: %v", err)
	}
	debug.DbgSlotStd("cos(2pi x/fT) =", out)
	debug.DbgChain("Chain before EvalCos :", eval, out)
	debug.PrintPrecStd("Prec out: ", want, out)

	maxErr := maxAbsErrReal(t, ecd, dec, out, want)
	t.Logf("EvalCos: fT=%g max err=%.3e", fT, maxErr)
	if maxErr > 1e-3 {
		t.Fatalf("EvalCos error %.3e too large", maxErr)
	}
}

func TestEvalSin(t *testing.T) {
	params, ecd, enc, dec, eval := newTrigoContext(t)
	ct, x := encryptInput(t, params, ecd, enc)

	fT := testPeriod * math.Pow(2, float64(testR))
	want := make([]float64, len(x))
	for s := range x {
		want[s] = math.Sin(2 * math.Pi * x[s] / fT)
	}

	debug.DbgSlotStd("input x =", ct)
	debug.DbgChain("Chain before EvalSin :", eval, ct)
	debug.PrintPrecStd("Prec in: ", x, ct)
	out, err := EvalSin(params, eval, ct, testK, testPeriod, testR, testDegree)
	if err != nil {
		t.Fatalf("EvalSin: %v", err)
	}
	debug.DbgSlotStd("sin(2pi x/fT) =", out)
	debug.DbgChain("Chain before EvalSin :", eval, out)
	debug.PrintPrecStd("Prec out: ", want, out)

	maxErr := maxAbsErrReal(t, ecd, dec, out, want)
	t.Logf("EvalSin: fT=%g max err=%.3e", fT, maxErr)
	if maxErr > 1e-3 {
		t.Fatalf("EvalSin error %.3e too large", maxErr)
	}
}

func TestEvalExpDirect(t *testing.T) {
	params, ecd, enc, dec, eval := newTrigoContext(t)
	ct, x := encryptInput(t, params, ecd, enc)

	// after r squarings the base frequency fT = period*2^r reaches 1/period.
	want := make([]complex128, len(x))
	for s := range x {
		want[s] = cmplx.Exp(complex(0, 2*math.Pi*x[s]/testPeriod))
	}

	debug.DbgSlotStd("input x =", ct)
	debug.DbgChain("Chain before EvalExpDirect :", eval, ct)
	debug.PrintPrecStd("Prec in: ", x, ct)
	out, err := EvalExpDirect(params, eval, ct, testK, testPeriod, testR, testDegree)
	if err != nil {
		t.Fatalf("EvalExpDirect: %v", err)
	}
	debug.DbgSlotStd("exp(2pi i x/T) =", out)
	debug.DbgChain("Chain after EvalExpDirect :", eval, out)
	debug.PrintPrecStd("Prec out: ", want, out)

	got := make([]complex128, out.Slots())
	if err := ecd.Decode(dec.DecryptNew(out), got); err != nil {
		t.Fatalf("decode: %v", err)
	}
	maxErr := 0.0
	for s := range want {
		if e := cmplx.Abs(got[s] - want[s]); e > maxErr {
			maxErr = e
		}
	}
	t.Logf("EvalExpDirect: max err=%.3e", maxErr)
	if maxErr > 1e-3 {
		t.Fatalf("EvalExpDirect error %.3e too large", maxErr)
	}
}

// maxAbsErrReal decrypts ct and returns the max |Re(got[s]) - want[s]|.
func maxAbsErrReal(t *testing.T, ecd *ckks.Encoder, dec *rlwe.Decryptor, ct *rlwe.Ciphertext, want []float64) float64 {
	t.Helper()
	got := make([]complex128, ct.Slots())
	if err := ecd.Decode(dec.DecryptNew(ct), got); err != nil {
		t.Fatalf("decode: %v", err)
	}
	maxErr := 0.0
	for s := range want {
		if e := math.Abs(real(got[s]) - want[s]); e > maxErr {
			maxErr = e
		}
	}
	return maxErr
}
