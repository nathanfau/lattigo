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

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/convctx"
	"github.com/tuneinsight/lattigo/v6/nathanfau/trigo"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const (
	// evalCosDeg is the Chebyshev degree of EvalCos and EvalSin.
	evalCosDeg = 31
	// evalExpR is the number of double-angle squarings applied after the trig step.
	evalExpR = 3
)

// step1 performs step 1 of Algorithm 1, that is CTS after ModRaise after STC(ct),
// and recombines the real and imaginary halves into one complex ciphertext.
func step1(eval *bootstrapping.Evaluator, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	//debug.DbgChain("  step1 in       :", eval.Evaluator, ct)

	ctSTC, err := eval.SlotsToCoeffs(ct, nil)
	if err != nil {
		return nil, fmt.Errorf("step1 SlotsToCoeffs: %w", err)
	}
	//debug.DbgCoeff("  after STC      :", ctSTC)
	//debug.DbgChain("  after STC      :", eval.Evaluator, ctSTC)

	ctSD, _, err := eval.ScaleDown(ctSTC)
	if err != nil {
		return nil, fmt.Errorf("step1 ScaleDown: %w", err)
	}
	//debug.DbgCoeff("  after ScaleDown:", ctSD)
	//debug.DbgChain("  after ScaleDown:", eval.Evaluator, ctSD)

	ctMU, err := eval.ModUp(ctSD)
	if err != nil {
		return nil, fmt.Errorf("step1 ModUp: %w", err)
	}
	//debug.DbgChain("  after ModUp    :", eval.Evaluator, ctMU)

	ctReal, ctImag, err := eval.CoeffsToSlots(ctMU)
	if err != nil {
		return nil, fmt.Errorf("after CoeffsToSlots: %w", err)
	}
	if ctImag == nil {
		return nil, fmt.Errorf("step1 CoeffsToSlots returned nil ctImag (sparse CTS), recombination impossible")
	}
	//debug.DbgSlotStd("  after CTS real :", ctReal)
	//debug.DbgSlotStd("  after CTS imag :", ctImag)
	//debug.DbgChain("  after CTS      :", eval.Evaluator, ctReal)

	e := eval.Evaluator
	ib := ctImag.CopyNew()
	if err := e.Mul(ib, complex(0, 1), ib); err != nil {
		return nil, fmt.Errorf("step1 Mul i: %w", err)
	}

	ct1 := ctReal.CopyNew()
	utils.AlignLevels(e, ct1, ib)
	if err := e.Add(ct1, ib, ct1); err != nil {
		return nil, fmt.Errorf("step1 Add: %w", err)
	}
	//debug.DbgSlotStd("  step1 out (ct1):", ct1)
	//debug.DbgChain("  step1 out (ct1):", eval.Evaluator, ct1)
	return ct1, nil
}

// Algo1 implements the 12 lines of Algorithm 1. Given a ciphertext packing integers
// m_s in {0, ..., t-1} (t = 2^k), it returns the real and imaginary parts of the
// t-th roots of unity exp(2*pi*i*m_s/t).
func Algo1(eval *bootstrapping.Evaluator, sw *convctx.CtxSwitcher, ct *rlwe.Ciphertext, k int) (ctreal, ctimag *rlwe.Ciphertext, err error) {
	// 1. ct1 <- CTS after ModRaise after STC(ct), recombined into one complex ct.
	//fmt.Println("---- Algo1 line 1: step1 (STC, ModRaise, CTS) ----")
	ct1, err := step1(eval, ct)
	if err != nil {
		return nil, nil, fmt.Errorf("line 1: %w", err)
	}
	//debug.DbgSlotStd("ct1 (Std)      :", ct1)
	//debug.DbgChain("ct1 (Std)      :", eval.Evaluator, ct1)

	// 2. ct2 <- Conv complex-to-real(ct1).
	//fmt.Println("---- Algo1 line 2: StandardToCI ----")
	ct2, err := sw.StandardToCI(ct1)
	if err != nil {
		return nil, nil, fmt.Errorf("line 2 StandardToCI: %w", err)
	}
	//debug.DbgSlotCI("ct2 (CI)       :", ct2)
	//debug.DbgChain("ct2 (CI)       :", sw.EvalCI, ct2)

	// 3-4. ctcos <- EvalCos(ct2), ctsin <- EvalSin(ct2), at base frequency (no squarings).
	t := 1 << k
	period := 1.0 / float64(t)
	//fmt.Println("---- Algo1 lines 3-4: EvalCos / EvalSin ----")
	ctcos, err := trigo.EvalCos(sw.CiP, sw.EvalCI, ct2, 1, period, evalExpR, evalCosDeg)
	if err != nil {
		return nil, nil, fmt.Errorf("line 3 EvalCos: %w", err)
	}
	//debug.DbgSlotCI("ctcos (CI)     :", ctcos)
	//debug.DbgChain("ctcos (CI)     :", sw.EvalCI, ctcos)
	ctsin, err := trigo.EvalSin(sw.CiP, sw.EvalCI, ct2, 1, period, evalExpR, evalCosDeg)
	if err != nil {
		return nil, nil, fmt.Errorf("line 4 EvalSin: %w", err)
	}
	//debug.DbgSlotCI("ctsin (CI)     :", ctsin)
	//debug.DbgChain("ctsin (CI)     :", sw.EvalCI, ctsin)

	// 5-6. Back to the complex (standard) context.
	//fmt.Println("---- Algo1 lines 5-6: CIToStandard ----")
	ctcosC, err := sw.CIToStandard(ctcos)
	if err != nil {
		return nil, nil, fmt.Errorf("line 5 CIToStandard cos: %w", err)
	}
	//debug.DbgSlotStd("ctcosC (Std)   :", ctcosC)
	//debug.DbgChain("ctcosC (Std)   :", eval.Evaluator, ctcosC)
	ctsinC, err := sw.CIToStandard(ctsin)
	if err != nil {
		return nil, nil, fmt.Errorf("line 6 CIToStandard sin: %w", err)
	}
	//debug.DbgSlotStd("ctsinC (Std)   :", ctsinC)
	//debug.DbgChain("ctsinC (Std)   :", eval.Evaluator, ctsinC)

	// 7-12. Re/Im extraction and recombination (extractExp is exactly these 6 lines).
	//fmt.Println("---- Algo1 lines 7-12: extractExp ----")
	ctreal, ctimag, err = extractExp(eval.Evaluator, ctcosC, ctsinC)
	if err != nil {
		return nil, nil, fmt.Errorf("lines 7-12 extractExp: %w", err)
	}
	//debug.DbgSlotStd("ctreal (Std)   :", ctreal)
	//debug.DbgSlotStd("ctimag (Std)   :", ctimag)
	//debug.DbgChain("ctreal (Std)   :", eval.Evaluator, ctreal)

	// Target frequency reached by r double-angle squarings (x2^r on the angle).
	//fmt.Println("---- Algo1 squareExp ----")
	if err = squareExp(eval.Evaluator, ctreal, ctimag, evalExpR); err != nil {
		return nil, nil, fmt.Errorf("squareExp: %w", err)
	}
	//debug.DbgSlotStd("ctreal sq (Std):", ctreal)
	//debug.DbgSlotStd("ctimag sq (Std):", ctimag)
	//debug.DbgChain("ctreal sq (Std):", eval.Evaluator, ctreal)
	return ctreal, ctimag, nil
}

// extractExp splits cos and sin into
//
//	ctReal = Re(cos) + i*Re(sin)
//	ctImag = Im(cos) + i*Im(sin)
func extractExp(e *ckks.Evaluator, ctCosC, ctSinC *rlwe.Ciphertext) (ctReal, ctImag *rlwe.Ciphertext, err error) {
	ctCosReal, err := ctRePart(e, ctCosC)
	if err != nil {
		return nil, nil, fmt.Errorf("Re cos: %w", err)
	}
	ctCosImag, err := ctImPart(e, ctCosC)
	if err != nil {
		return nil, nil, fmt.Errorf("Im cos: %w", err)
	}
	ctSinReal, err := ctRePart(e, ctSinC)
	if err != nil {
		return nil, nil, fmt.Errorf("Re sin: %w", err)
	}
	ctSinImag, err := ctImPart(e, ctSinC)
	if err != nil {
		return nil, nil, fmt.Errorf("Im sin: %w", err)
	}

	if ctReal, err = utils.CombineReIm(e, ctCosReal, ctSinReal); err != nil {
		return nil, nil, fmt.Errorf("combine ct_real: %w", err)
	}
	if ctImag, err = utils.CombineReIm(e, ctCosImag, ctSinImag); err != nil {
		return nil, nil, fmt.Errorf("combine ct_imag: %w", err)
	}
	return ctReal, ctImag, nil
}

// ctRePart returns Re(z) = (z + conj(z))/2.
func ctRePart(eval *ckks.Evaluator, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	conj := ct.CopyNew()
	if err := eval.Conjugate(conj, conj); err != nil {
		return nil, fmt.Errorf("conj: %w", err)
	}
	out := ct.CopyNew()
	if err := eval.Add(out, conj, out); err != nil {
		return nil, fmt.Errorf("add: %w", err)
	}
	if err := eval.Mul(out, complex(0.5, 0), out); err != nil {
		return nil, fmt.Errorf("mul 1/2: %w", err)
	}
	if err := eval.Rescale(out, out); err != nil {
		return nil, fmt.Errorf("rescale: %w", err)
	}
	return out, nil
}

// ctImPart returns Im(z) = (conj(z) - z)*i/2.
func ctImPart(eval *ckks.Evaluator, ct *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	conj := ct.CopyNew()
	if err := eval.Conjugate(conj, conj); err != nil {
		return nil, fmt.Errorf("conj: %w", err)
	}
	out := conj.CopyNew()
	if err := eval.Sub(out, ct, out); err != nil {
		return nil, fmt.Errorf("sub: %w", err)
	}
	if err := eval.Mul(out, complex(0, 0.5), out); err != nil {
		return nil, fmt.Errorf("mul i/2: %w", err)
	}
	if err := eval.Rescale(out, out); err != nil {
		return nil, fmt.Errorf("rescale: %w", err)
	}
	return out, nil
}

// squareExp applies r double-angle squarings to the real and imaginary parts,
// each squared independently (MulRelin then Rescale), to reach the target frequency.
func squareExp(e *ckks.Evaluator, ctreal, ctimag *rlwe.Ciphertext, r int) error {
	for i := 0; i < r; i++ {
		if err := e.MulRelin(ctreal, ctreal, ctreal); err != nil {
			return fmt.Errorf("squaring %d ctreal: %w", i+1, err)
		}
		if err := e.Rescale(ctreal, ctreal); err != nil {
			return fmt.Errorf("rescale squaring %d ctreal: %w", i+1, err)
		}
		if err := e.MulRelin(ctimag, ctimag, ctimag); err != nil {
			return fmt.Errorf("squaring %d ctimag: %w", i+1, err)
		}
		if err := e.Rescale(ctimag, ctimag); err != nil {
			return fmt.Errorf("rescale squaring %d ctimag: %w", i+1, err)
		}
	}
	return nil
}
