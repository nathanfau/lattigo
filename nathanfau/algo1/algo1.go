// Package algo1 implements Algorithm 1 of [BCKK25], the conjugate-invariant
// IntRootBoot (a bootstrap variant built on EvalCos, EvalSin and a CI domain
// switch).
//
// We name it "algo1" on purpose. Calling it "bbbts" would be misleading. For us
// BBBTS is the name of the (different) algorithm in [BKSS24], migrated in the
// bbbts package, so we keep the paper's "Algorithm 1" label here to avoid any
// confusion.
package algo1

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/convctx"
	"github.com/tuneinsight/lattigo/v6/nathanfau/trigo"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const (
	// evalCosDeg is the Chebyshev degree of EvalCos and EvalSin.
	evalCosDeg = 30
	// evalExpR is the number of double-angle squarings applied after the trig step.
	evalExpR = 3
	// entryTol is how close to 1, in bits, ExtractTo wants what ScaleDown leaves uncorrected. It
	// only has to separate an input brought to EntryScale (off by the scale's own rounding, far
	// below 2^-100) from one that was not (off by at least a prime's distance to a power of two).
	entryTol = 30
)

// step1 performs step 1 of Algorithm 1, that is CTS after ModRaise after STC(ct),
// and recombines the real and imaginary halves into one complex ciphertext. It also returns what
// ScaleDown could not correct (1 when ct arrives at EntryScale).
func step1(eval *bootstrapping.Evaluator, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, rlwe.Scale, error) {

	//debug.DbgSlotStd("step1 in:", ct)
	//debug.DbgChain("step1 in:", eval.Evaluator, ct)

	ctSTC, err := eval.SlotsToCoeffs(ct, nil)
	if err != nil {
		return nil, rlwe.Scale{}, fmt.Errorf("step1 SlotsToCoeffs: %w", err)
	}

	//debug.DbgCoeff("after STC:", ctSTC)
	//debug.DbgChain("after STC:", eval.Evaluator, ctSTC)

	ctSD, errScale, err := eval.ScaleDown(ctSTC)
	if err != nil {
		return nil, rlwe.Scale{}, fmt.Errorf("step1 ScaleDown: %w", err)
	}

	//debug.DbgCoeff("after ScaleDown:", ctSD)
	//debug.DbgChain("after ScaleDown:", eval.Evaluator, ctSD)

	ctMU, err := eval.ModUp(ctSD)
	if err != nil {
		return nil, rlwe.Scale{}, fmt.Errorf("step1 ModUp: %w", err)
	}

	//debug.DbgChain("after ModUp:", eval.Evaluator, ctMU)

	ctReal, ctImag, err := eval.CoeffsToSlots(ctMU)
	if err != nil {
		return nil, rlwe.Scale{}, fmt.Errorf("after CoeffsToSlots: %w", err)
	}
	if ctImag == nil {
		return nil, rlwe.Scale{}, fmt.Errorf("step1 CoeffsToSlots returned nil ctImag (sparse CTS), recombination impossible")
	}

	//debug.DbgSlotStd("after CTS real:", ctReal)
	//debug.DbgSlotStd("after CTS imag:", ctImag)
	//debug.DbgChain("after CTS:", eval.Evaluator, ctReal)

	ct1, err := utils.CombineReIm(eval.Evaluator, ctReal, ctImag)
	if err != nil {
		return nil, rlwe.Scale{}, fmt.Errorf("step1 CombineReIm: %w", err)
	}

	//debug.DbgSlotStd("step1 out (ct1):", ct1)
	//debug.DbgChain("step1 out (ct1):", eval.Evaluator, ct1)

	return ct1, *errScale, nil
}

// EntryScale is the scale Algo1 wants its input at: Q[0]/MessageRatio, the one ScaleDown brings
// the ciphertext to before ModUp. SlotsToCoeffs does not move the scale, so a ciphertext arriving
// at EntryScale goes through ScaleDown untouched. Any other scale is corrected by ScaleDown only up
// to an integer factor; the rest (its errScale) reaches the mod-1 as a relative error on the
// message, which Algo1, unlike EvalMod, has no way to divide out afterwards.
func EntryScale(eval *bootstrapping.Evaluator) rlwe.Scale {
	q0 := rlwe.NewScale(eval.BootstrappingParameters.Q()[0])
	return q0.Div(rlwe.NewScale(eval.Mod1Parameters.MessageRatio()))
}

// squareInput is the scale Resume's squarings must start from, at level lvl, to end exactly on
// out. A squaring maps s to s^2/q, q what its Rescale divides by, so over evalExpR = r of them
//
//	out = s^(2^r) / prod_{j<r} q_{lvl-j}^(2^(r-1-j))
//
// and s is the 2^r-th root of out times that product.
func squareInput(params ckks.Parameters, lvl int, out rlwe.Scale) rlwe.Scale {
	acc := new(big.Float).SetPrec(2 * rlwe.ScalePrecision).Set(&out.Value)
	for j := 0; j < evalExpR; j++ {
		q := utils.RescalingFactor(params, lvl-j)
		for e := 0; e < evalExpR-1-j; e++ {
			q = q.Mul(q)
		}
		acc.Mul(acc, &q.Value)
	}
	for j := 0; j < evalExpR; j++ {
		acc.Sqrt(acc)
	}
	return rlwe.NewScale(acc)
}

// Run implements the 12 lines of Algorithm 1. Given a ciphertext packing integers
// m_s in {0, ..., t-1} (t = 2^k), it returns the real and imaginary parts of the
// t-th roots of unity exp(2*pi*i*m_s/t).
func Run(eval *bootstrapping.Evaluator, sw *convctx.CtxSwitcher, ct *rlwe.Ciphertext, k int) (ctreal, ctimag *rlwe.Ciphertext, err error) {
	if ctreal, ctimag, err = Extract(eval, sw, ct, k); err != nil {
		return nil, nil, err
	}
	if err = Resume(eval, ctreal, ctimag); err != nil {
		return nil, nil, err
	}
	return ctreal, ctimag, nil
}

// Extract runs Algorithm 1 up to and including the Re/Im extraction (extractExp): the
// point where a packet's real part (m_A) and imaginary part (m_B) become two separate
// ciphertexts
func Extract(eval *bootstrapping.Evaluator, sw *convctx.CtxSwitcher, ct *rlwe.Ciphertext, k int) (ctreal, ctimag *rlwe.Ciphertext, err error) {
	return ExtractTo(eval, sw, ct, k, rlwe.Scale{})
}

// ExtractTo is Extract for a ciphertext that comes from, and goes back to, a zone of smaller primes
// (params2.Shape.LogZone), whose scale out is not the bootstrap's. Its two edges:
//
//   - in: ct must arrive at EntryScale (convctx.CIToStandardTo takes it there). ExtractTo fails if
//     ScaleDown still had something to correct, since Algo1 could not undo it.
//   - out: extractExp multiplies by a constant, so it can land on any scale for free. It lands on
//     the one from which Resume's squarings end exactly on out, where the bit extraction runs.
//
// A zero out is Extract: no check, and extractExp keeps the scale it is given.
func ExtractTo(eval *bootstrapping.Evaluator, sw *convctx.CtxSwitcher, ct *rlwe.Ciphertext, k int, out rlwe.Scale) (ctreal, ctimag *rlwe.Ciphertext, err error) {
	zone := out.Value.Sign() != 0

	// 1. ct1 <- CTS after ModRaise after STC(ct), recombined into one complex ct

	//fmt.Println("---- Algo1 line 1: step1 (STC, ModRaise, CTS) ----")

	ct1, errScale, err := step1(eval, ct)
	if err != nil {
		return nil, nil, fmt.Errorf("line 1: %w", err)
	}
	if off := errScale.Log2Delta(rlwe.NewScale(1)); zone && off < entryTol {
		return nil, nil, fmt.Errorf("line 1: ScaleDown left an error of 2^-%.1f on the message, want below 2^-%d: the input is at 2^%.6f, not at EntryScale 2^%.6f",
			off, entryTol, ct.Scale.Log2(), EntryScale(eval).Log2())
	}

	//debug.DbgSlotStd("ct1 (Std):", ct1)
	//debug.DbgChain("ct1 (Std):", eval.Evaluator, ct1)

	// 2. ct2 <- Conv_{Cplx->Real} (ct1)

	//fmt.Println("---- Algo1 line 2: StandardToCI ----")

	ct2, err := sw.StandardToCI(ct1)
	if err != nil {
		return nil, nil, fmt.Errorf("line 2 StandardToCI: %w", err)
	}

	//debug.DbgSlotCI("ct2 (CI):", ct2)
	//debug.DbgChain("ct2 (CI):", sw.EvalCI, ct2)

	// 3-4. ctcos <- EvalCos(ct2), ctsin <- EvalSin(ct2), at base frequency (no squarings).

	//fmt.Println("---- Algo1 lines 3-4: EvalCos / EvalSin ----")

	t := 1 << k
	period := 1.0 / float64(t)

	ctcos, err := trigo.EvalCos(sw.CiP, sw.EvalCI, ct2, 1, period, evalExpR, evalCosDeg)
	if err != nil {
		return nil, nil, fmt.Errorf("line 3 EvalCos: %w", err)
	}

	//debug.DbgSlotCI("ctcos (CI):", ctcos)
	//debug.DbgChain("ctcos (CI):", sw.EvalCI, ctcos)

	ctsin, err := trigo.EvalSin(sw.CiP, sw.EvalCI, ct2, 1, period, evalExpR, evalCosDeg)
	if err != nil {
		return nil, nil, fmt.Errorf("line 4 EvalSin: %w", err)
	}

	//debug.DbgSlotCI("ctsin (CI):", ctsin)
	//debug.DbgChain("ctsin (CI):", sw.EvalCI, ctsin)

	// 5-6. Back to the complex (std) context

	//fmt.Println("---- Algo1 lines 5-6: CIToStandard ----")

	ctcosC, err := sw.CIToStandard(ctcos)
	if err != nil {
		return nil, nil, fmt.Errorf("line 5 CIToStandard cos: %w", err)
	}

	//debug.DbgSlotStd("ctcosC (Std):", ctcosC)
	//debug.DbgChain("ctcosC (Std):", eval.Evaluator, ctcosC)

	ctsinC, err := sw.CIToStandard(ctsin)
	if err != nil {
		return nil, nil, fmt.Errorf("line 6 CIToStandard sin: %w", err)
	}

	//debug.DbgSlotStd("ctsinC (Std):", ctsinC)
	//debug.DbgChain("ctsinC (Std):", eval.Evaluator, ctsinC)

	// 7-12. Re/Im extraction  ---- Algo1 pauses here in the refresh. ----

	//fmt.Println("---- Algo1 lines 7-12: extractExp ----")

	var mid rlwe.Scale // zero: extractExp keeps the scale it is given
	if zone {
		// extractExp rescales once, so Resume starts one level below its inputs.
		mid = squareInput(eval.BootstrappingParameters, min(ctcosC.Level(), ctsinC.Level())-1, out)
	}
	ctreal, ctimag, err = extractExp(eval.Evaluator, ctcosC, ctsinC, mid)
	if err != nil {
		return nil, nil, fmt.Errorf("lines 7-12 extractExp: %w", err)
	}

	//debug.DbgSlotStd("ctreal (Std):", ctreal)
	//debug.DbgSlotStd("ctimag (Std):", ctimag)
	//debug.DbgChain("ctreal/ctimag chain (Std):", eval.Evaluator, ctreal)

	return ctreal, ctimag, nil
}

// Resume finishes Algorithm 1 after the pause: the evalExpR double-angle squarings that
// bring the extracted Re/Im parts from the base frequency to the target t-th roots of unity.
func Resume(eval *bootstrapping.Evaluator, ctreal, ctimag *rlwe.Ciphertext) error {
	if err := squareExp(eval.Evaluator, ctreal, ctimag, evalExpR); err != nil {
		return fmt.Errorf("squareExp: %w", err)
	}

	//debug.DbgSlotStd("ctreal squarred (Std):", ctreal)
	//debug.DbgSlotStd("ctimag squarred (Std):", ctimag)
	//debug.DbgChain("ctreal squarred (Std):", eval.Evaluator, ctreal)

	return nil
}

// extractExp splits cos and sin into
//
//	ctReal = Re(cos) + i*Re(sin) = (C + conj(C) + i(S + conj(S))) / 2
//	ctImag = Im(cos) + i*Im(sin) = (i(conj(C) - C) - (conj(S) - S)) / 2
//
// Both come out at the scale out, or at the scale they are given when out is zero: the / 2 is the
// multiplication by a constant that sets it.
func extractExp(e *ckks.Evaluator, ctCosC, ctSinC *rlwe.Ciphertext, out rlwe.Scale) (ctReal, ctImag *rlwe.Ciphertext, err error) {
	half := func(ct *rlwe.Ciphertext) error {
		if out.Value.Sign() != 0 {
			return utils.MulRescaleTo(e, ct, 0.5, out)
		}
		if err := e.Mul(ct, complex(0.5, 0), ct); err != nil {
			return fmt.Errorf("*1/2: %w", err)
		}
		if err := e.Rescale(ct, ct); err != nil {
			return fmt.Errorf("rescale: %w", err)
		}
		return nil
	}

	Cc := ctCosC.CopyNew()
	if err = e.Conjugate(Cc, Cc); err != nil {
		return nil, nil, fmt.Errorf("conj cos: %w", err)
	}
	Sc := ctSinC.CopyNew()
	if err = e.Conjugate(Sc, Sc); err != nil {
		return nil, nil, fmt.Errorf("conj sin: %w", err)
	}

	// ctReal = (C + conj(C) + i(S + conj(S))) / 2
	reC := ctCosC.CopyNew()
	if err = e.Add(reC, Cc, reC); err != nil { // C + conj(C)
		return nil, nil, fmt.Errorf("real C+conj(C): %w", err)
	}
	reS := ctSinC.CopyNew()
	if err = e.Add(reS, Sc, reS); err != nil { // S + conj(S)
		return nil, nil, fmt.Errorf("real S+conj(S): %w", err)
	}
	if err = e.Mul(reS, complex(0, 1), reS); err != nil { // i(S + conj(S))
		return nil, nil, fmt.Errorf("real *i: %w", err)
	}
	utils.AlignLevels(e, reC, reS)
	if err = e.Add(reC, reS, reC); err != nil {
		return nil, nil, fmt.Errorf("real sum: %w", err)
	}
	if err = half(reC); err != nil { // / 2
		return nil, nil, fmt.Errorf("real: %w", err)
	}
	ctReal = reC

	// ctImag = (i(conj(C) - C) - (conj(S) - S)) / 2
	imC := Cc.CopyNew()
	if err = e.Sub(imC, ctCosC, imC); err != nil { // conj(C) - C
		return nil, nil, fmt.Errorf("imag conj(C)-C: %w", err)
	}
	if err = e.Mul(imC, complex(0, 1), imC); err != nil { // i(conj(C) - C)
		return nil, nil, fmt.Errorf("imag *i: %w", err)
	}
	imS := Sc.CopyNew()
	if err = e.Sub(imS, ctSinC, imS); err != nil { // conj(S) - S
		return nil, nil, fmt.Errorf("imag conj(S)-S: %w", err)
	}
	utils.AlignLevels(e, imC, imS)
	if err = e.Sub(imC, imS, imC); err != nil { // i(conj(C)-C) - (conj(S)-S)
		return nil, nil, fmt.Errorf("imag diff: %w", err)
	}
	if err = half(imC); err != nil { // / 2
		return nil, nil, fmt.Errorf("imag: %w", err)
	}
	ctImag = imC

	return ctReal, ctImag, nil
}

// square applies one double-angle squaring in place (MulRelin then Rescale).
func square(e *ckks.Evaluator, ct *rlwe.Ciphertext) error {
	if err := e.MulRelin(ct, ct, ct); err != nil {
		return fmt.Errorf("MulRelin: %w", err)
	}
	if err := e.Rescale(ct, ct); err != nil {
		return fmt.Errorf("Rescale: %w", err)
	}
	return nil
}

// squareExp applies r double-angle squarings to the real and imaginary parts,
// each squared independently (MulRelin then Rescale), to reach the target frequency.
func squareExp(e *ckks.Evaluator, ctreal, ctimag *rlwe.Ciphertext, r int) error {
	for i := 0; i < r; i++ {
		if err := square(e, ctreal); err != nil {
			return fmt.Errorf("squaring %d ctreal: %w", i+1, err)
		}
		if err := square(e, ctimag); err != nil {
			return fmt.Errorf("squaring %d ctimag: %w", i+1, err)
		}
	}
	return nil
}
