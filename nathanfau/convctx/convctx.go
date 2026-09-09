// Package convctx converts ciphertexts between the conjugate-invariant (CI) and the
// standard (Std) CKKS domains at full capacity (BCKK25, Fig. 1 and 2).
//
//	CI  (deg N, N real slots)      : [x_0 .. x_{N/2-1}, y_0 .. y_{N/2-1}]
//	Std (deg N, N/2 complex slots) : [z_0 .. z_{N/2-1}]   with z_j = x_j + i*y_j
//
// The first half of the CI real slots hold the real parts, the second half the imaginary parts
package convctx

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// CtxSwitcher holds all the material (keys, masks, evaluators) for the CI <-> Std conversion.
type CtxSwitcher struct {
	CiP       ckks.Parameters
	StdBigP   ckks.Parameters
	StdSmallP ckks.Parameters

	SkCI    *rlwe.SecretKey
	SkBig   *rlwe.SecretKey // iota(SkSmall), not an independent key
	SkSmall *rlwe.SecretKey

	nHalf int

	EvalCI *ckks.Evaluator
	eval   *ckks.Evaluator

	sw         ckks.DomainSwitcher
	evkStdToCI *rlwe.EvaluationKey // iota(s) -> s~; the other way round, sw.RealToComplex holds it

	stdRingQ *ring.Ring // deg 2N, standard
	ciRingQ  *ring.Ring // deg N, conjugate invariant

	autRot  []uint64 // sigma_{2N+1}, the N/2-slot rotation, as a coefficient permutation
	autConj []uint64 // sigma_{-1}, the conjugation, as a coefficient permutation

	// One mask per level, encoded at that level's rescaling factor so a conversion does not move
	// the scale. Built on demand by maskAt.
	ecdBig  *ckks.Encoder
	ptMask1 map[int]*rlwe.Plaintext // [1.., i..]
	ptMask2 map[int]*rlwe.Plaintext // [1.., -i..]
}

// NewCtxSwitcher builds the machinery around a Std deg-N ring (typically the bootstrapping
// ring) and its secret key skSmall, so that ciphertexts of that ring (e.g. the output of
// CoeffsToSlots) are directly convertible. The CI deg-N ring and the Std deg-2N working ring
// are derived from the same literal (same moduli Q, P); the CI secret is generated here.
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

	// SkBig = iota(SkSmall). MapSmallDimensionToLargerDimensionNTT is lattigo's iota in the NTT
	// domain -- GenEvaluationKey uses it to embed the smaller of two keys -- and the P part is
	// rebuilt from Q the way GenEvaluationKeysForRingSwapNew does.
	c.SkBig = rlwe.NewSecretKey(c.StdBigP)
	ring.MapSmallDimensionToLargerDimensionNTT(skSmall.Value.Q, c.SkBig.Value.Q)
	if c.StdBigP.PCount() != 0 {
		buffQ := c.StdBigP.RingQ().NewPoly()
		rlwe.ExtendBasisSmallNormAndCenterNTTMontgomery(c.StdBigP.RingQ(), c.StdBigP.RingP(), c.SkBig.Value.Q, buffQ, c.SkBig.Value.P)
	}

	c.SkCI = rlwe.NewKeyGenerator(c.CiP).GenSecretKeyNew()

	// The only two key-switching keys of the package: s~ -> iota(s), which the DomainSwitcher uses
	// in RealToComplex, and iota(s) -> s~, which StandardToCI applies itself because it has to
	// place it before the masking.
	evkC2R, evkR2C := kgenBig.GenEvaluationKeysForRingSwapNew(c.SkBig, c.SkCI)
	c.evkStdToCI = evkC2R
	if c.sw, err = ckks.NewDomainSwitcher(c.StdBigP, evkC2R, evkR2C); err != nil {
		return nil, fmt.Errorf("NewDomainSwitcher: %w", err)
	}

	// The two automorphisms, as coefficient permutations rather than as key switches.
	// GaloisElement(N/2) = 5^{N/2} mod 4N = 2N+1, the sigma of Fig. 1.
	c.stdRingQ = c.StdBigP.RingQ()
	if c.ciRingQ, err = c.stdRingQ.ConjugateInvariantRing(); err != nil {
		return nil, fmt.Errorf("ConjugateInvariantRing: %w", err)
	}
	if c.autRot, err = ring.AutomorphismNTTIndex(c.stdRingQ.N(), c.stdRingQ.NthRoot(), c.StdBigP.GaloisElement(c.nHalf)); err != nil {
		return nil, fmt.Errorf("AutomorphismNTTIndex rot: %w", err)
	}
	if c.autConj, err = ring.AutomorphismNTTIndex(c.stdRingQ.N(), c.stdRingQ.NthRoot(), c.stdRingQ.NthRoot()-1); err != nil {
		return nil, fmt.Errorf("AutomorphismNTTIndex conj: %w", err)
	}

	c.eval = ckks.NewEvaluator(c.StdBigP, nil)

	rlkCI := rlwe.NewKeyGenerator(c.CiP).GenRelinearizationKeyNew(c.SkCI)
	c.EvalCI = ckks.NewEvaluator(c.CiP, rlwe.NewMemEvaluationKeySet(rlkCI))

	c.ecdBig = ckks.NewEncoder(c.StdBigP)
	c.ptMask1 = map[int]*rlwe.Plaintext{}
	c.ptMask2 = map[int]*rlwe.Plaintext{}

	return c, nil
}

// maskAt returns the mask of Fig. 1 step (2), or of Fig. 2 step (3) when second, for a ciphertext
// at level lvl.
func (c *CtxSwitcher) maskAt(lvl int, second bool) (*rlwe.Plaintext, error) {

	cache := c.ptMask1
	if second {
		cache = c.ptMask2
	}
	if pt, ok := cache[lvl]; ok {
		return pt, nil
	}

	N := c.StdBigP.MaxSlots()
	half := N / 2
	m := make([]complex128, N)
	for i := 0; i < half; i++ {
		m[i] = 1
	}
	v := 1i
	if second {
		v = -1i
	}
	for i := half; i < N; i++ {
		m[i] = v
	}

	// At what Rescale will divide by, so the conversion does not move the scale.
	scale := utils.RescalingFactor(c.StdBigP, lvl)
	if second {
		scale = scale.Div(rlwe.NewScale(2))
	}

	pt := ckks.NewPlaintext(c.StdBigP, lvl)
	pt.Scale = scale
	if err := c.ecdBig.Encode(m, pt); err != nil {
		return nil, fmt.Errorf("encode mask (second=%t) at level %d: %w", second, lvl, err)
	}
	cache[lvl] = pt
	return pt, nil
}

// CIToStandard is Fig. 1
func (c *CtxSwitcher) CIToStandard(ctCI *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	if ctCI.Level() < 1 {
		return nil, fmt.Errorf("CIToStandard: input at level %d, the masking needs a prime to rescale", ctCI.Level())
	}

	ctBig := ckks.NewCiphertext(c.StdBigP, 1, ctCI.Level())
	if err := c.sw.RealToComplex(c.eval, ctCI, ctBig); err != nil {
		return nil, fmt.Errorf("RealToComplex: %w", err)
	}

	pt, err := c.maskAt(ctBig.Level(), false)
	if err != nil {
		return nil, err
	}
	if err := c.eval.Mul(ctBig, pt, ctBig); err != nil {
		return nil, fmt.Errorf("mul mask1: %w", err)
	}
	if err := c.eval.Rescale(ctBig, ctBig); err != nil {
		return nil, fmt.Errorf("rescale mask1: %w", err)
	}

	level := ctBig.Level()
	ringQ := c.stdRingQ.AtLevel(level)
	ctRot := ckks.NewCiphertext(c.StdBigP, 1, level)
	for i := range ctBig.Value {
		ringQ.AutomorphismNTTWithIndex(ctBig.Value[i], c.autRot, ctRot.Value[i])
		ringQ.Add(ctBig.Value[i], ctRot.Value[i], ctBig.Value[i])
	}

	ctSmall := ckks.NewCiphertext(c.StdSmallP, 1, level)
	rlwe.SwitchCiphertextRingDegreeNTT(ctBig.El(), ringQ, ctSmall.El())
	ctSmall.LogDimensions = c.StdSmallP.LogMaxDimensions()
	return ctSmall, nil
}

// StandardToCI is Fig. 2
func (c *CtxSwitcher) StandardToCI(ctStd *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	if ctStd.Level() < 1 {
		return nil, fmt.Errorf("StandardToCI: input at level %d, the masking needs a prime to rescale", ctStd.Level())
	}

	ctBig := ckks.NewCiphertext(c.StdBigP, 1, ctStd.Level())
	rlwe.SwitchCiphertextRingDegreeNTT(ctStd.El(), nil, ctBig.El())
	ctBig.LogDimensions = c.StdBigP.LogMaxDimensions()

	if err := c.eval.ApplyEvaluationKey(ctBig, c.evkStdToCI, ctBig); err != nil {
		return nil, fmt.Errorf("key switch iota(s) -> s~: %w", err)
	}

	pt, err := c.maskAt(ctBig.Level(), true)
	if err != nil {
		return nil, err
	}
	if err := c.eval.Mul(ctBig, pt, ctBig); err != nil {
		return nil, fmt.Errorf("mul mask2: %w", err)
	}
	if err := c.eval.Rescale(ctBig, ctBig); err != nil {
		return nil, fmt.Errorf("rescale mask2: %w", err)
	}

	level := ctBig.Level()
	ciRingQ := c.ciRingQ.AtLevel(level)
	ctCI := ckks.NewCiphertext(c.CiP, 1, level)
	for i := range ctBig.Value {
		ciRingQ.FoldStandardToConjugateInvariant(ctBig.Value[i], c.autConj, ctCI.Value[i])
	}
	*ctCI.MetaData = *ctBig.MetaData
	ctCI.LogDimensions = c.CiP.LogMaxDimensions()

	// FLAG SCALE. The fold doubled the encoded value, so the scale doubles with it.
	ctCI.Scale = ctBig.Scale.Mul(rlwe.NewScale(2))
	return ctCI, nil
}
