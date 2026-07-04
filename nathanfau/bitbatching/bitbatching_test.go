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
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const testK = 4

// newTestContext builds a standard (complex) CKKS context with relinearization and
// complex-conjugation keys, deep enough to run BitExtract for testK bits.
func newTestContext(t *testing.T) (ckks.Parameters, *ckks.Encoder, *rlwe.Encryptor, *rlwe.Decryptor, *ckks.Evaluator) {
	t.Helper()
	logQ := []int{55}
	for i := 0; i < 6; i++ {
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

// encryptBit encrypts a full slot vector
func encryptBit(t *testing.T, params ckks.Parameters, ecd *ckks.Encoder, enc *rlwe.Encryptor, vals []float64) *rlwe.Ciphertext {
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

// decodeReal decrypts ct and returns the real part of each slot.
func decodeReal(t *testing.T, ecd *ckks.Encoder, dec *rlwe.Decryptor, ct *rlwe.Ciphertext) []float64 {
	t.Helper()
	out := make([]complex128, ct.Slots())
	if err := ecd.Decode(dec.DecryptNew(ct), out); err != nil {
		t.Fatalf("decode: %v", err)
	}
	re := make([]float64, len(out))
	for i, v := range out {
		re[i] = real(v)
	}
	return re
}

// TestBitPack: encrypt the k bit-planes of random integers m, pack them, and check
// the packed ciphertext decrypts to m = sum_i b_i * 2^i.
func TestBitPack(t *testing.T) {
	params, ecd, enc, dec, eval := newTestContext(t)
	k := testK
	n := params.MaxSlots()
	rng := rand.New(rand.NewSource(time.Now().UnixNano())) // *random* seed for testing

	m := make([]int, n)
	planes := make([][]float64, k)
	for i := range planes {
		planes[i] = make([]float64, n)
	}
	for s := 0; s < n; s++ {
		m[s] = rng.Intn(1 << k)
		for i := 0; i < k; i++ {
			planes[i][s] = float64((m[s] >> i) & 1)
		}
	}

	for i := 0; i < k; i++ {
		fmt.Printf("planes[%d] = %f\n", i, planes[i][:2])
	}

	bits := make([]*rlwe.Ciphertext, k)
	for i := 0; i < k; i++ {
		bits[i] = encryptBit(t, params, ecd, enc, planes[i])
	}

	for i := 0; i < k; i++ {
		debug.DbgSlotStd("encryptedplanes = ", bits[i])
	}

	ctPack, err := BitPack(eval, bits)
	if err != nil {
		t.Fatalf("bitPack: %v", err)
	}

	debug.DbgSlotStd("BitPack =", ctPack)

	got := decodeReal(t, ecd, dec, ctPack)
	maxErr := 0.0
	for s := 0; s < n; s++ {
		if e := math.Abs(got[s] - float64(m[s])); e > maxErr {
			maxErr = e
		}
	}
	t.Logf("bitPack max error = %.3e", maxErr)
	if maxErr > 0.1 {
		t.Fatalf("bitPack: max error %.3e too large", maxErr)
	}
}

// TestBitExtract: encode omega^m directly and check if BitExtract recovers each bit of m.
func TestBitExtract(t *testing.T) {
	params, ecd, enc, dec, eval := newTestContext(t)
	k := testK
	n := params.MaxSlots()
	rng := rand.New(rand.NewSource(time.Now().UnixNano())) // *random* seed for testing

	m := make([]int, n)
	vals := make([]complex128, n)
	base := 2 * math.Pi / float64(int(1)<<k)
	for s := 0; s < n; s++ {
		m[s] = rng.Intn(1 << k)
		vals[s] = cmplx.Exp(complex(0, base*float64(m[s])))
	}

	fmt.Printf("vals = %f\n", vals[:2])

	pt := ckks.NewPlaintext(params, params.MaxLevel())
	if err := ecd.Encode(vals, pt); err != nil {
		t.Fatalf("encode omega^m: %v", err)
	}
	ctExp, err := enc.EncryptNew(pt)
	if err != nil {
		t.Fatalf("encrypt omega^m: %v", err)
	}

	debug.DbgSlotStd("omega^m =", ctExp)
	debug.DbgChain("Chain before extraction :", eval, ctExp)

	bits, err := BitExtract(params, eval, ctExp, k)
	if err != nil {
		t.Fatalf("BitExtract: %v", err)
	}
	if len(bits) != k {
		t.Fatalf("BitExtract returned %d bits, want %d", len(bits), k)
	}

	for i := 0; i < 4; i++ {
		debug.DbgSlotStd("bitsExtracted = ", bits[i])
	}

	for i := 0; i < 4; i++ {
		debug.DbgChain("Chain after extraction :", eval, bits[i])
	}

	maxErr := 0.0
	for i := 0; i < k; i++ {
		got := decodeReal(t, ecd, dec, bits[i])
		for s := 0; s < n; s++ {
			want := float64((m[s] >> i) & 1)
			if e := math.Abs(got[s] - want); e > maxErr {
				maxErr = e
			}
			if math.Round(got[s]) != want {
				t.Fatalf("bit %d slot %d: got %.4f, want %.0f (m=%d)", i, s, got[s], want, m[s])
			}
		}
	}
	t.Logf("BitExtract max bit error = %.3e", maxErr)
}
