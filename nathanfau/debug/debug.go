package debug

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

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

// dbgPrefix prints the aligned "label [ctx, logN= ,  lv= , sc= ]" prefix shared by DbgPrinters
// so their labels and outputs line up in the same columns.
func dbgPrefix(label, ctx string, ct *rlwe.Ciphertext) {
	fmt.Printf("  %-30s  [%-3s, logN=%-2d, lv=%-2d sc=2^%-9.3f]   ", label, ctx, ct.LogN(), ct.Level(), math.Log2(ct.Scale.Float64()))
}

// dbgSlot decrypts ct with the given decoding context and prints its first slots under label.
// The three exported variants below differ only by that context.
func dbgSlot(label, ctx string, ecd *ckks.Encoder, dec *rlwe.Decryptor, ct *rlwe.Ciphertext) {
	if ecd == nil || dec == nil {
		return
	}
	buf := make([]complex128, 2)
	if err := ecd.Decode(dec.DecryptNew(ct), buf); err != nil {
		dbgPrefix(label, ctx, ct)
		fmt.Printf("decode error: %v\n", err)
		return
	}
	dbgPrefix(label, ctx, ct)
	for _, v := range buf {
		fmt.Printf("%s ", fmtVal(v, false))
	}
	fmt.Println("...")
}

// DbgSlotStd decrypts ct and prints its first slots in the Std context.
func DbgSlotStd(label string, ct *rlwe.Ciphertext) {
	dbgSlot(label, "Std", EncStd, DecStd, ct)
}

// DbgSlotCI is the conjugate-invariant (CI) variant of DbgSlotStd.
func DbgSlotCI(label string, ct *rlwe.Ciphertext) {
	dbgSlot(label, "CI", EncCI, DecCI, ct)
}

// DbgSlotBig is the Big-context variant of DbgSlotStd.
func DbgSlotBig(label string, ct *rlwe.Ciphertext) {
	dbgSlot(label, "Big", EncBig, DecBig, ct)
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

// printPrec prints the precision statistics of ct against want, in the given decoding context.
// realOnly keeps the REAL column only, for a conjugate-invariant (real) ciphertext.
func printPrec(label string, want any, ct *rlwe.Ciphertext, params ckks.Parameters, ecd *ckks.Encoder, dec *rlwe.Decryptor, realOnly bool) {
	if ecd == nil || dec == nil {
		return
	}
	prec := ckks.GetPrecisionStats(params, ecd, dec, want, ct, 0, false)
	table := prec.String()
	if realOnly {
		table = precStringReal(prec)
	}
	fmt.Printf(" %s \n%s", label, table)
}

// PrintPrecStd prints the precision statistics of ct against the wanted values (Std context).
// want may be []float64 or []complex128 (anything GetPrecisionStats accepts).
func PrintPrecStd(label string, want interface{}, ct *rlwe.Ciphertext) {
	printPrec(label, want, ct, ParamsStd, EncStd, DecStd, false)
}

// DbgBitStd traces one bit-ciphertext of a bit-sliced circuit in the Std context under a single
// label: its first slots, its modulus chain, then its precision against want.
func DbgBitStd(label string, eval *ckks.Evaluator, ct *rlwe.Ciphertext, want interface{}) {
	DbgSlotStd("ct_"+label+" =", ct)
	DbgChain("chain ct_"+label+" :", eval, ct)
	PrintPrecStd("prec ct_"+label+" :", want, ct)
}

// PrintPrecCI is the conjugate-invariant (CI) variant of PrintPrecStd.
func PrintPrecCI(label string, want interface{}, ct *rlwe.Ciphertext) {
	printPrec(label, want, ct, ParamsCI, EncCI, DecCI, true)
}

// precStringReal formats a PrecisionStats table with the REAL column only.
func precStringReal(p ckks.PrecisionStats) string {
	f := func(v float64) string { return fmt.Sprintf("%.2f", v) }
	return "\n" + utils.BoxTable("Log2",
		[]utils.TableGroup{{Name: "REAL", Cols: []string{"Prec", "Err"}}},
		[]string{"MIN", "MAX", "AVG", "MED", "STD"},
		[][]string{
			{f(p.MINLog2Prec.Real), f(p.MINLog2Err.Real)},
			{f(p.MAXLog2Prec.Real), f(p.MAXLog2Err.Real)},
			{f(p.AVGLog2Prec.Real), f(p.AVGLog2Err.Real)},
			{f(p.MEDLog2Prec.Real), f(p.MEDLog2Err.Real)},
			{f(p.STDLog2Prec.Real), f(p.STDLog2Err.Real)},
		})
}

// PrecCt pairs a ciphertext with the values it should decrypt to (want is []float64 or
// []complex128). Name is an optional label for the worst-slot line (falls back to "ct <index>").
type PrecCt struct {
	Name string
	Want any
	Ct   *rlwe.Ciphertext
}

// PrecPoolStd prints the SAME statistics table as PrintPrecStd
// but POOLED over several ciphertexts
func PrecPoolStd(label string, entries ...PrecCt) {
	precPool(label, entries, ParamsStd, EncStd, DecStd, false)
}

// PrecPoolCI is the conjugate-invariant (real) variant of PrecPoolStd. It prints only the REAL
// column: a CI ciphertext decodes to real values, so IMAG/L2 precision is meaningless.
func PrecPoolCI(label string, entries ...PrecCt) {
	precPool(label, entries, ParamsCI, EncCI, DecCI, true)
}

func precPool(label string, entries []PrecCt, params ckks.Parameters, ecd *ckks.Encoder, dec *rlwe.Decryptor, realOnly bool) {
	if ecd == nil || dec == nil || len(entries) == 0 {
		return
	}
	var wantAll, haveAll []complex128
	stats := make([]utils.BitStats, len(entries))
	decoded := make([][]complex128, len(entries))
	for c, e := range entries {
		buf := make([]complex128, e.Ct.Slots())
		if err := ecd.Decode(dec.DecryptNew(e.Ct), buf); err != nil {
			fmt.Printf(" %s : decode error (ct %d): %v\n", label, c, err)
			return
		}
		wc := utils.WantToCplx(e.Want, len(buf))
		stats[c] = utils.BitDistanceOf(buf, e.Want, 0)
		decoded[c] = buf
		wantAll = append(wantAll, wc...)
		haveAll = append(haveAll, buf[:len(wc)]...)
	}
	agg, worstCt := utils.WorstOf(stats)
	bestCt := 0
	for c, st := range stats {
		if st.Best == agg.Best {
			bestCt = c
			break
		}
	}
	worstSlot, bestSlot := agg.WorstSlot, agg.BestSlot
	worstErr, bestErr := agg.Worst, agg.Best
	worstGot, worstWant := decoded[worstCt][worstSlot], utils.WantToCplx(entries[worstCt].Want, len(decoded[worstCt]))[worstSlot]
	bestGot, bestWant := decoded[bestCt][bestSlot], utils.WantToCplx(entries[bestCt].Want, len(decoded[bestCt]))[bestSlot]
	prec := ckks.GetPrecisionStats(params, ecd, dec, wantAll, haveAll, 0, false)
	table := prec.String()
	if realOnly {
		table = precStringReal(prec)
	}
	fmt.Printf(" %s (%d ct, %d slots) \n%s", label, len(entries), len(wantAll), table)
	fmt.Printf(" WORST: %s slot %d  err=%s, margin to flip = %.4f  got=%s want=%s\n",
		ctName(entries, worstCt), worstSlot, fmtErr(worstErr), 0.5-agg.Worst, fmtVal(worstGot, realOnly), fmtVal(worstWant, realOnly))
	fmt.Printf(" BEST : %s slot %d  err=%s  got=%s want=%s\n",
		ctName(entries, bestCt), bestSlot, fmtErr(bestErr), fmtVal(bestGot, realOnly), fmtVal(bestWant, realOnly))
}

// ctName is the label of entry c (falls back to "ct <index>").
func ctName(entries []PrecCt, c int) string {
	if entries[c].Name != "" {
		return entries[c].Name
	}
	return fmt.Sprintf("ct %d", c)
}

// fmtErr formats an absolute error as "e (2^log2)", or "0 (exact)" when it is exactly zero.
func fmtErr(er float64) string {
	if er <= 0 {
		return "0 (exact)"
	}
	return fmt.Sprintf("%.3e (2^%.2f)", er, math.Log2(er))
}

// fmtVal prints a slot value: just the real part for a CI (real) pool, else the full complex value.
func fmtVal(v complex128, realOnly bool) string {
	if realOnly {
		return fmt.Sprintf("%+.15f", real(v))
	}
	return fmt.Sprintf("(%+.15f%+.15fi)", real(v), imag(v))
}

// DbgChain prints the remaining chain of primes of the ct,
func DbgChain(label string, eval *ckks.Evaluator, ct *rlwe.Ciphertext) {
	Q := eval.GetRLWEParameters().Q()
	logs := make([]int, ct.Level()+1)
	for i := range logs {
		logs[i] = int(math.Round(math.Log2(float64(Q[i]))))
	}
	dbgPrefix(label, "Chn", ct)
	fmt.Printf("%d\n", logs)
}

// DbgParams prints a summary of parameters
func DbgParams(name string, p ckks.Parameters) {
	fmt.Printf("  RingType        : %s\n", p.RingType())
	fmt.Printf("  LogN            : %d  (ring degree 2^%d)\n", p.LogN(), p.LogN())
	fmt.Printf("  LogSlots        : %d  (MaxSlots %d)\n", p.LogMaxSlots(), p.MaxSlots())
	fmt.Printf("  LogDefaultScale : %d\n", p.LogDefaultScale())
	fmt.Printf("  Levels          : %d primes (MaxLevel %d)\n", p.QCount(), p.MaxLevel())
	fmt.Printf("  Xs HammingWeight: %d\n", p.XsHammingWeight())
	fmt.Printf("  logQ  (per q_i) : %d  = %.1f\n", p.LogQi(), p.LogQ())
	fmt.Printf("  logP  (per p_i) : %d  = %.1f\n", p.LogPi(), p.LogP())
	fmt.Printf("  logQP           : %.1f\n", p.LogQP())
}
