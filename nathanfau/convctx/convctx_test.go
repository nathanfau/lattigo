package convctx

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// TestCtxSwitchRoundTrip encrypts a complex vector z in the Std deg-N ring, converts it to
// the CI ring (Re/Im split into the two halves of the real slots), and back, checking both
// the split and the full round trip.
func TestCtxSwitchRoundTrip(t *testing.T) {
	// Build a CI ring first so the primes are NTT-friendly ( = 1 mod (4N))
	ciBase, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            12,
		LogQ:            []int{55, 45, 45, 45},
		LogP:            []int{55},
		LogDefaultScale: 45,
		RingType:        ring.ConjugateInvariant,
	})
	if err != nil {
		t.Fatalf("ci params: %v", err)
	}
	stdLit := ciBase.ParametersLiteral()
	stdLit.RingType = ring.Standard
	stdSmall, err := ckks.NewParametersFromLiteral(stdLit)
	if err != nil {
		t.Fatalf("std params: %v", err)
	}

	skSmall := rlwe.NewKeyGenerator(stdSmall).GenSecretKeyNew()

	c, err := NewCtxSwitcher(stdSmall, skSmall)
	if err != nil {
		t.Fatalf("NewCtxSwitcher: %v", err)
	}

	L := stdSmall.MaxSlots()
	half := c.CiP.MaxSlots() / 2
	if half != L {
		t.Fatalf("layout mismatch: CI half=%d, Std slots=%d", half, L)
	}

	ecdStd := ckks.NewEncoder(stdSmall)
	encStd := rlwe.NewEncryptor(stdSmall, skSmall)
	decStd := rlwe.NewDecryptor(stdSmall, skSmall)
	evalStd := ckks.NewEvaluator(stdSmall, nil)

	ecdCI := ckks.NewEncoder(c.CiP)
	decCI := rlwe.NewDecryptor(c.CiP, c.SkCI)

	debug.EncStd, debug.DecStd, debug.ParamsStd = ecdStd, decStd, stdSmall
	debug.EncCI, debug.DecCI, debug.ParamsCI = ecdCI, decCI, c.CiP

	// Random complex input z
	rng := rand.New(rand.NewSource(time.Now().UnixNano())) // *random* seed for testing
	z := make([]complex128, L)
	for i := range z {
		z[i] = complex(rng.Float64()*2-1, rng.Float64()*2-1)
	}

	fmt.Printf("cplx values: %f ...\n", z[:2])

	// Expected CI layout after StandardToCI
	wantCI := make([]float64, c.CiP.MaxSlots())
	for j := 0; j < L; j++ {
		wantCI[j] = real(z[j])
		wantCI[half+j] = imag(z[j])
	}

	pt := ckks.NewPlaintext(stdSmall, stdSmall.MaxLevel())
	if err := ecdStd.Encode(z, pt); err != nil {
		t.Fatalf("encode: %v", err)
	}
	ctStd, err := encStd.EncryptNew(pt)
	if err != nil {
		t.Fatalf("encrypt: %v", err)
	}
	debug.DbgSlotStd("ct_z =", ctStd)
	debug.DbgChain("chain ct_z :", evalStd, ctStd)
	debug.PrintPrecStd("prec ct_z :", z, ctStd)

	// Std -> CI
	ctCI, err := c.StandardToCI(ctStd)
	if err != nil {
		t.Fatalf("StandardToCI: %v", err)
	}
	debug.DbgSlotCI("ct_CI =", ctCI)
	fmt.Printf("(only real values are seen above...)\n")
	debug.DbgChain("chain ct_CI:", c.EvalCI, ctCI)
	debug.PrintPrecCI("prec ct_CI :", wantCI, ctCI)

	ci := make([]float64, c.CiP.MaxSlots())
	if err := ecdCI.Decode(decCI.DecryptNew(ctCI), ci); err != nil {
		t.Fatalf("decode CI: %v", err)
	}
	maxErrSplit := 0.0
	for j := 0; j < L; j++ {
		if e := math.Abs(ci[j] - real(z[j])); e > maxErrSplit {
			maxErrSplit = e
		}
		if e := math.Abs(ci[half+j] - imag(z[j])); e > maxErrSplit {
			maxErrSplit = e
		}
	}
	t.Logf("Std->CI split max err = %.3e", maxErrSplit)
	if maxErrSplit > 1e-2 {
		t.Fatalf("Std->CI split error %.3e too large", maxErrSplit)
	}

	// CI -> Std
	ctStd2, err := c.CIToStandard(ctCI)
	if err != nil {
		t.Fatalf("CIToStandard: %v", err)
	}
	debug.DbgSlotStd("ct_z' =", ctStd2)
	debug.DbgChain("chain ct_z' :", evalStd, ctStd2)
	debug.PrintPrecStd("prec ct_z' :", z, ctStd2)

	got := make([]complex128, ctStd2.Slots())
	if err := ecdStd.Decode(decStd.DecryptNew(ctStd2), got); err != nil {
		t.Fatalf("decode Std: %v", err)
	}
	maxErrRT := 0.0
	for j := 0; j < L; j++ {
		if e := math.Abs(real(got[j])-real(z[j])) + math.Abs(imag(got[j])-imag(z[j])); e > maxErrRT {
			maxErrRT = e
		}
	}
	t.Logf("Std->CI->Std round-trip max err = %.3e", maxErrRT)
	if maxErrRT > 1e-2 {
		t.Fatalf("round-trip error %.3e too large", maxErrRT)
	}
}
