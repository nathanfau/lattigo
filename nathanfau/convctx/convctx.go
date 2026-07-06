// Package convctx converts ciphertexts between the conjugate-invariant (CI) and the
// standard (Std) CKKS domains at full capacity (BCKK25, Fig. 1 and 2).
//
//	CI  (deg N, N real slots)      : [x_0 .. x_{N/2-1}, y_0 .. y_{N/2-1}]
//	Std (deg N, N/2 complex slots) : [z_0 .. z_{N/2-1}]   with z_j = x_j + i*y_j
//
// The first half of the CI real slots hold the real parts, the second half the imaginary
// parts. An intermediate Std ring of degree 2N (N complex slots) bridges the two:
// RealToComplex/ComplexToReal cross CI <-> Std-deg-2N, and ApplyEvaluationKey changes the
// degree 2N <-> N.
package convctx

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// CtxSwitcher holds all the material (keys, masks, evaluators) for the CI <-> Std conversion.
type CtxSwitcher struct {
	CiP       ckks.Parameters
	StdBigP   ckks.Parameters
	StdSmallP ckks.Parameters

	SkCI    *rlwe.SecretKey
	SkBig   *rlwe.SecretKey
	SkSmall *rlwe.SecretKey

	nHalf int

	EvalCI *ckks.Evaluator
	eval   *ckks.Evaluator

	sw ckks.DomainSwitcher

	evkBigToSmall *rlwe.EvaluationKey // 2N -> N projection
	evkSmallToBig *rlwe.EvaluationKey // N -> 2N embedding

	ptMask1 *rlwe.Plaintext // [1.., i..]
	ptMask2 *rlwe.Plaintext // [1/2.., -i/2..]
}

// NewCtxSwitcher builds the machinery around a Std deg-N ring (typically the bootstrapping
// ring) and its secret key skSmall, so that ciphertexts of that ring (e.g. the output of
// CoeffsToSlots) are directly convertible. The CI deg-N ring and the Std deg-2N working ring
// are derived from the same literal (same moduli Q, P); their secrets are generated here.
func NewCtxSwitcher(stdSmallP ckks.Parameters, skSmall *rlwe.SecretKey) (*CtxSwitcher, error) {

	c := &CtxSwitcher{StdSmallP: stdSmallP, SkSmall: skSmall}

	ciLit := stdSmallP.ParametersLiteral()
	ciLit.RingType = ring.ConjugateInvariant
	var err error
	if c.CiP, err = ckks.NewParametersFromLiteral(ciLit); err != nil {
		return nil, fmt.Errorf("CiP: %w", err)
	}

	bigLit := stdSmallP.ParametersLiteral()
	bigLit.LogN = stdSmallP.LogN() + 1
	bigLit.RingType = ring.Standard
	if c.StdBigP, err = ckks.NewParametersFromLiteral(bigLit); err != nil {
		return nil, fmt.Errorf("StdBigP: %w", err)
	}

	c.nHalf = c.StdBigP.MaxSlots() / 2

	kgenBig := rlwe.NewKeyGenerator(c.StdBigP)
	c.SkBig = kgenBig.GenSecretKeyNew()
	c.SkCI = rlwe.NewKeyGenerator(c.CiP).GenSecretKeyNew()

	evkC2R, evkR2C := kgenBig.GenEvaluationKeysForRingSwapNew(c.SkBig, c.SkCI)
	if c.sw, err = ckks.NewDomainSwitcher(c.StdBigP, evkC2R, evkR2C); err != nil {
		return nil, fmt.Errorf("NewDomainSwitcher: %w", err)
	}

	c.evkBigToSmall = kgenBig.GenEvaluationKeyNew(c.SkBig, c.SkSmall)
	c.evkSmallToBig = kgenBig.GenEvaluationKeyNew(c.SkSmall, c.SkBig)

	galRot := kgenBig.GenGaloisKeyNew(c.StdBigP.GaloisElement(c.nHalf), c.SkBig)
	galConj := kgenBig.GenGaloisKeyNew(c.StdBigP.GaloisElementOrderTwoOrthogonalSubgroup(), c.SkBig)
	evkSet := rlwe.NewMemEvaluationKeySet(nil, galRot, galConj)

	c.eval = ckks.NewEvaluator(c.StdBigP, evkSet)

	rlkCI := rlwe.NewKeyGenerator(c.CiP).GenRelinearizationKeyNew(c.SkCI)
	c.EvalCI = ckks.NewEvaluator(c.CiP, rlwe.NewMemEvaluationKeySet(rlkCI))

	N := c.StdBigP.MaxSlots()
	half := N / 2
	ecdBig := ckks.NewEncoder(c.StdBigP)

	m1 := make([]complex128, N)
	for i := 0; i < half; i++ {
		m1[i] = 1
	}
	for i := half; i < N; i++ {
		m1[i] = 1i
	}
	c.ptMask1 = ckks.NewPlaintext(c.StdBigP, c.StdBigP.MaxLevel())
	if err = ecdBig.Encode(m1, c.ptMask1); err != nil {
		return nil, fmt.Errorf("encode mask1: %w", err)
	}

	m2 := make([]complex128, N)
	for i := 0; i < half; i++ {
		m2[i] = 0.5
	}
	for i := half; i < N; i++ {
		m2[i] = -0.5i
	}
	c.ptMask2 = ckks.NewPlaintext(c.StdBigP, c.StdBigP.MaxLevel())

	/*


		FLAG


	*/
	// Compensate the folding factor 2: ComplexToReal (FoldStandardToConjugateInvariant)
	// cost = ~ 1 bit of precision
	// need to analyse this further
	c.ptMask2.Scale = c.ptMask2.Scale.Div(rlwe.NewScale(2))
	if err = ecdBig.Encode(m2, c.ptMask2); err != nil {
		return nil, fmt.Errorf("encode mask2: %w", err)
	}

	return c, nil
}

// CIToStandard is Fig. 1
func (c *CtxSwitcher) CIToStandard(ctCI *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	ctBig := ckks.NewCiphertext(c.StdBigP, 1, ctCI.Level())
	if err := c.sw.RealToComplex(c.eval, ctCI, ctBig); err != nil {
		return nil, fmt.Errorf("RealToComplex: %w", err)
	}

	if err := c.eval.Mul(ctBig, c.ptMask1, ctBig); err != nil {
		return nil, fmt.Errorf("mul mask1: %w", err)
	}
	if err := c.eval.Rescale(ctBig, ctBig); err != nil {
		return nil, fmt.Errorf("rescale mask1: %w", err)
	}

	ctRot := ckks.NewCiphertext(c.StdBigP, 1, ctBig.Level())
	if err := c.eval.Rotate(ctBig, c.nHalf, ctRot); err != nil {
		return nil, fmt.Errorf("rotate: %w", err)
	}
	if err := c.eval.Add(ctBig, ctRot, ctBig); err != nil {
		return nil, fmt.Errorf("add: %w", err)
	}

	ctSmall := ckks.NewCiphertext(c.StdSmallP, 1, ctBig.Level())
	if err := c.eval.ApplyEvaluationKey(ctBig, c.evkBigToSmall, ctSmall); err != nil {
		return nil, fmt.Errorf("project 2N->N: %w", err)
	}

	ctSmall.LogDimensions = c.StdSmallP.LogMaxDimensions()
	return ctSmall, nil
}

// StandardToCI is Fig. 2
func (c *CtxSwitcher) StandardToCI(ctStd *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	ctBig := ckks.NewCiphertext(c.StdBigP, 1, ctStd.Level())
	if err := c.eval.ApplyEvaluationKey(ctStd, c.evkSmallToBig, ctBig); err != nil {
		return nil, fmt.Errorf("embed N->2N: %w", err)
	}
	ctBig.LogDimensions = c.StdBigP.LogMaxDimensions()

	ctMasked := ckks.NewCiphertext(c.StdBigP, 1, ctBig.Level())
	if err := c.eval.Mul(ctBig, c.ptMask2, ctMasked); err != nil {
		return nil, fmt.Errorf("mul mask2: %w", err)
	}
	if err := c.eval.Rescale(ctMasked, ctMasked); err != nil {
		return nil, fmt.Errorf("rescale mask2: %w", err)
	}

	ctConj := ckks.NewCiphertext(c.StdBigP, 1, ctMasked.Level())
	if err := c.eval.Conjugate(ctMasked, ctConj); err != nil {
		return nil, fmt.Errorf("conjugate: %w", err)
	}
	if err := c.eval.Add(ctMasked, ctConj, ctMasked); err != nil {
		return nil, fmt.Errorf("add: %w", err)
	}

	ctCI := ckks.NewCiphertext(c.CiP, 1, ctMasked.Level())
	if err := c.sw.ComplexToReal(c.eval, ctMasked, ctCI); err != nil {
		return nil, fmt.Errorf("ComplexToReal: %w", err)
	}

	return ctCI, nil
}
