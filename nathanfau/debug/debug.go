package debug

import (
	"fmt"
	"math"
	"strings"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// TimeSince prints the elapsed time since t0 under the given label + returns it
func TimeSince(label string, t0 time.Time) time.Duration {
	d := time.Since(t0)
	fmt.Printf("%-20s %s\n", label, d.Round(time.Microsecond))
	return d
}

// Decoding context (encoder, decryptor, parameters).
// Std is for the Standard context (Complex context)
// CI if for Conjugate-Invariant context (Real context)
// Big is for Big Standard context (LogN*2 Complex context)
var EncStd *ckks.Encoder
var DecStd *rlwe.Decryptor
var ParamsStd ckks.Parameters

var EncCI *ckks.Encoder
var DecCI *rlwe.Decryptor
var ParamsCI ckks.Parameters

var EncBig *ckks.Encoder
var DecBig *rlwe.Decryptor

// Slot decrypts ct and prints its first slots in the Std context.
func Slot(label string, ct *rlwe.Ciphertext) {
	if DecStd == nil || EncStd == nil {
		return
	}
	n := 2
	buf := make([]complex128, n)
	if err := EncStd.Decode(DecStd.DecryptNew(ct), buf); err != nil {
		fmt.Printf("  %-30s  error: %v\n", label, err)
		return
	}
	fmt.Printf("  %-30s  [Std logN=%d slots=%d]  lv=%-2d  sc=2^%.1f   ", label, ct.LogN(), ct.Slots(), ct.Level(), math.Log2(ct.Scale.Float64()))
	for _, v := range buf[:n] {
		fmt.Printf("(%+.15f%+.15fi) ", real(v), imag(v))
	}
	fmt.Println("...")
}

// SlotCI is the conjugate-invariant (CI) variant of Slot.
func SlotCI(label string, ct *rlwe.Ciphertext) {
	if DecCI == nil || EncCI == nil {
		return
	}
	n := 2
	buf := make([]complex128, n)
	if err := EncCI.Decode(DecCI.DecryptNew(ct), buf); err != nil {
		fmt.Printf("  %-30s  error: %v\n", label, err)
		return
	}
	fmt.Printf("  %-30s  [CI  logN=%d slots=%d]  lv=%-2d  sc=2^%.1f   ", label, ct.LogN(), ct.Slots(), ct.Level(), math.Log2(ct.Scale.Float64()))
	for _, v := range buf[:n] {
		fmt.Printf("(%+.15f%+.15fi) ", real(v), imag(v))
	}
	fmt.Println("...")
}

// SlotBig is the Big-context variant of Slot.
func SlotBig(label string, ct *rlwe.Ciphertext) {
	if DecBig == nil || EncBig == nil {
		return
	}
	n := 2
	buf := make([]complex128, n)
	if err := EncBig.Decode(DecBig.DecryptNew(ct), buf); err != nil {
		fmt.Printf("  %-30s  error: %v\n", label, err)
		return
	}
	fmt.Printf("  %-30s  [Std logN=%d slots=%d]  lv=%-2d  sc=2^%.1f   ", label, ct.LogN(), ct.Slots(), ct.Level(), math.Log2(ct.Scale.Float64()))
	for _, v := range buf[:n] {
		fmt.Printf("(%+.15f%+.15fi) ", real(v), imag(v))
	}
	fmt.Println("...")
}

// Coeff prints the first centered coefficients of a ciphertext's decrypted plaintext.
func Coeff(label string, ct *rlwe.Ciphertext) {
	if DecStd == nil {
		return
	}
	const n = 2
	pt := DecStd.DecryptNew(ct)
	q0 := ParamsStd.Q()[0]
	half := q0 >> 1
	scale := ct.Scale.Float64()
	fmt.Printf("  %-30s  lv=%-2d  sc=2^%.1f   coeffs: ", label, ct.Level(), math.Log2(scale))
	for i := range n {
		c := pt.Value.Coeffs[0][i]
		var v int64
		if c > half {
			v = int64(c) - int64(q0)
		} else {
			v = int64(c)
		}
		fmt.Printf("%+.15f ", float64(v)/scale)
	}
	fmt.Println("...")
}

// PrintPrec prints the precision statistics of ct against the wanted values (Std context).
func PrintPrec(label string, want []float64, ct *rlwe.Ciphertext) {
	if DecStd == nil || EncStd == nil {
		return
	}
	prec := ckks.GetPrecisionStats(ParamsStd, EncStd, DecStd, want, ct, 0, false)
	fmt.Printf(" %s \n%s", label, prec.String())
}

// PrintPrecCI is the conjugate-invariant (CI) variant of PrintPrec.
func PrintPrecCI(label string, want []float64, ct *rlwe.Ciphertext) {
	if DecCI == nil || EncCI == nil {
		return
	}
	prec := ckks.GetPrecisionStats(ParamsCI, EncCI, DecCI, want, ct, 0, false)
	fmt.Printf(" %s \n%s", label, prec.String())
}

// Chain prints the remaining chain of primes of the ct,
func Chain(label string, eval *ckks.Evaluator, ct *rlwe.Ciphertext) {
	Q := eval.GetRLWEParameters().Q()
	var sb strings.Builder
	for i := 0; i <= ct.Level(); i++ {
		if i > 0 {
			sb.WriteByte(' ')
		}
		fmt.Fprintf(&sb, "%d", int(math.Round(math.Log2(float64(Q[i])))))
	}
	fmt.Printf("       %-34s lv=%-2d sc=2^%-5.1f [%s]\n",
		label, ct.Level(), math.Log2(ct.Scale.Float64()), sb.String())
}
