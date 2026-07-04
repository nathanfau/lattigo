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

// dbgPrefix prints the aligned "label [ctx, logN= ] lv= sc=" prefix shared by DbgPrinters
// so their labels and outputs line up in the same columns.
func dbgPrefix(label, ctx string, ct *rlwe.Ciphertext) {
	fmt.Printf("  %-30s  [%-3s, logN=%-2d, lv=%-2d sc=2^%-6.1f]   ", label, ctx, ct.LogN(), ct.Level(), math.Log2(ct.Scale.Float64()))
}

// DbgSlotStd decrypts ct and prints its first slots in the Std context.
func DbgSlotStd(label string, ct *rlwe.Ciphertext) {
	if DecStd == nil || EncStd == nil {
		return
	}
	n := 2
	buf := make([]complex128, n)
	if err := EncStd.Decode(DecStd.DecryptNew(ct), buf); err != nil {
		dbgPrefix(label, "Std", ct)
		fmt.Printf("decode error: %v\n", err)
		return
	}
	dbgPrefix(label, "Std", ct)
	for _, v := range buf[:n] {
		fmt.Printf("(%+.15f%+.15fi) ", real(v), imag(v))
	}
	fmt.Println("...")
}

// DbgSlotCI is the conjugate-invariant (CI) variant of DbgSlotStd.
func DbgSlotCI(label string, ct *rlwe.Ciphertext) {
	if DecCI == nil || EncCI == nil {
		return
	}
	n := 2
	buf := make([]complex128, n)
	if err := EncCI.Decode(DecCI.DecryptNew(ct), buf); err != nil {
		dbgPrefix(label, "CI", ct)
		fmt.Printf("decode error: %v\n", err)
		return
	}
	dbgPrefix(label, "CI", ct)
	for _, v := range buf[:n] {
		fmt.Printf("(%+.15f%+.15fi) ", real(v), imag(v))
	}
	fmt.Println("...")
}

// DbgSlotBig is the Big-context variant of DbgSlotStd.
func DbgSlotBig(label string, ct *rlwe.Ciphertext) {
	if DecBig == nil || EncBig == nil {
		return
	}
	n := 2
	buf := make([]complex128, n)
	if err := EncBig.Decode(DecBig.DecryptNew(ct), buf); err != nil {
		dbgPrefix(label, "Big", ct)
		fmt.Printf("decode error: %v\n", err)
		return
	}
	dbgPrefix(label, "Big", ct)
	for _, v := range buf[:n] {
		fmt.Printf("(%+.15f%+.15fi) ", real(v), imag(v))
	}
	fmt.Println("...")
}

// DbgCoeff prints the first coefficients of a ciphertext's decrypted plaintext.
func DbgCoeff(label string, ct *rlwe.Ciphertext) {
	if DecStd == nil {
		return
	}
	const n = 2
	pt := DecStd.DecryptNew(ct)
	q0 := ParamsStd.Q()[0]
	half := q0 >> 1
	scale := ct.Scale.Float64()
	dbgPrefix(label, "Cff", ct)
	fmt.Print("coeffs: ")
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

// PrintPrecStd prints the precision statistics of ct against the wanted values (Std context).
// want may be []float64 or []complex128 (anything GetPrecisionStats accepts).
func PrintPrecStd(label string, want interface{}, ct *rlwe.Ciphertext) {
	if DecStd == nil || EncStd == nil {
		return
	}
	prec := ckks.GetPrecisionStats(ParamsStd, EncStd, DecStd, want, ct, 0, false)
	fmt.Printf(" %s \n%s", label, prec.String())
}

// PrintPrecCI is the conjugate-invariant (CI) variant of PrintPrecStd.
func PrintPrecCI(label string, want interface{}, ct *rlwe.Ciphertext) {
	if DecCI == nil || EncCI == nil {
		return
	}
	prec := ckks.GetPrecisionStats(ParamsCI, EncCI, DecCI, want, ct, 0, false)
	fmt.Printf(" %s \n%s", label, prec.String())
}

// DbgChain prints the remaining chain of primes of the ct,
func DbgChain(label string, eval *ckks.Evaluator, ct *rlwe.Ciphertext) {
	Q := eval.GetRLWEParameters().Q()
	var sb strings.Builder
	for i := 0; i <= ct.Level(); i++ {
		if i > 0 {
			sb.WriteByte(' ')
		}
		fmt.Fprintf(&sb, "%d", int(math.Round(math.Log2(float64(Q[i])))))
	}
	dbgPrefix(label, "Chn", ct)
	fmt.Printf("[%s]\n", sb.String())
}
