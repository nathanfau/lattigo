// Package transciphering wires the Algo1 IntRootBoot refresh into a homomorphic AES
// middle-round pipeline so an encrypted AES state can be advanced round by round under CKKS.
//
// One "middle" round (neither the first AddRoundKey-only round nor the last MixColumns-less
// round) is:
//
//	SubBytes -> refresh(Algo1) -> ShiftRows -> MixColumns -> AddRoundKey -> Cleaning
//
// The AES circuit runs in the conjugate-invariant (CI, real) context of the CtxSwitcher; the
// refresh round-trips each bit CI -> Std, packs k bits per ciphertext, bootstraps with Algo1,
// extracts the bits back, and returns to CI. SubBytes and MixColumns versions are selectable.
package transciphering

import (
	"fmt"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
	"github.com/tuneinsight/lattigo/v6/nathanfau/algo1"
	"github.com/tuneinsight/lattigo/v6/nathanfau/bitbatching"
	"github.com/tuneinsight/lattigo/v6/nathanfau/cleaning"
	"github.com/tuneinsight/lattigo/v6/nathanfau/convctx"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"

	"math/big"
)

// Context bundles every piece needed to run the AES + refresh pipeline: the Std bootstrapping
// parameters/evaluator, the CtxSwitcher (CI <-> Std), the CI-context AES evaluator, and the CI
// encoder/encryptor/decryptor for the state.
type Context struct {
	Params ckks.Parameters // Std bootstrapping / residual parameters
	Eval   *bootstrapping.Evaluator
	Sw     *convctx.CtxSwitcher
	AE     *aes.Evaluator // AES round functions, CI context (sw.EvalCI)

	EcdCI *ckks.Encoder
	EncCI *rlwe.Encryptor
	DecCI *rlwe.Decryptor

	K      int        // bits packed per ciphertext (t = 2^K)
	Target int        // Algo1 input level (S2C LevelQ)
	Canon  rlwe.Scale // canonical output scale after refresh
}

// NewContext builds the whole machinery on the algo1 moduli chain (same setup as the algo1
// test). The AES state lives in the CI ring (sw.CiP); the bootstrap runs in the Std ring.
func NewContext(logN, k int, noTrick bool, mcVersion int) (*Context, error) {
	params, btpParams, err := buildBtpParams(logN, k, mcVersion)
	if err != nil {
		return nil, fmt.Errorf("buildBtpParams: %w", err)
	}

	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()

	fmt.Println("BTS KeyGen ...")
	t0 := time.Now()
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		return nil, fmt.Errorf("GenEvaluationKeys: %w", err)
	}
	eval, err := bootstrapping.NewEvaluator(btpParams, evk)
	if err != nil {
		return nil, fmt.Errorf("bootstrapping NewEvaluator: %w", err)
	}
	fmt.Printf("Done !(%s)\n", time.Since(t0).Round(time.Millisecond))

	fmt.Println("CtxSwitcher KeyGen ...")
	t0 = time.Now()
	sw, err := convctx.NewCtxSwitcher(params, sk)
	if err != nil {
		return nil, fmt.Errorf("NewCtxSwitcher: %w", err)
	}
	fmt.Printf("Done !(%s)\n", time.Since(t0).Round(time.Millisecond))

	ae := aes.NewEvaluator(sw.EvalCI)
	if noTrick {
		ae = aes.NewEvaluatorNoTrick(sw.EvalCI)
	}

	c := &Context{
		Params: params,
		Eval:   eval,
		Sw:     sw,
		AE:     ae,
		EcdCI:  ckks.NewEncoder(sw.CiP),
		EncCI:  rlwe.NewEncryptor(sw.CiP, sw.SkCI),
		DecCI:  rlwe.NewDecryptor(sw.CiP, sw.SkCI),
		K:      k,
		Target: eval.SlotsToCoeffsParameters.LevelQ,
		Canon:  sw.CiP.DefaultScale(),
	}

	// Debug decoding contexts: Std (residual) and CI (AES circuit).
	debug.EncStd, debug.DecStd, debug.ParamsStd = ckks.NewEncoder(params), rlwe.NewDecryptor(params, sk), params
	debug.EncCI, debug.DecCI, debug.ParamsCI = c.EcdCI, c.DecCI, sw.CiP

	return c, nil
}

// buildBtpParams builds the Std CKKS parameters and bootstrapping parameters on the algo1
// moduli chain. The chain is annotated per pipeline step; the user tunes the q_i.
func buildBtpParams(logN, k, mcVersion int) (ckks.Parameters, bootstrapping.Parameters, error) {
	logQ := []int{42}
	logQ = append(logQ, 60)         // SlotsToCoeffs
	logQ = append(logQ, 38)         // Conv_{Real->Cplx}
	logQ = append(logQ, 38, 38, 38) // SubBytes
	logQ = append(logQ, 38, 38)     // cleaning
	logQ = append(logQ, 38)         // AddRoundKey
	logQ = append(logQ, 38, 38, 38) // MixColumns (V2 depth 3)
	if mcVersion == 1 {
		logQ = append(logQ, 38) // MixColumnsV1 (depth 4) consumes 1 more level
	}
	// refresh
	logQ = append(logQ, 38)                 // Conv_{Cplx->Real}
	logQ = append(logQ, 38, 38, 38, 38, 38) // 5 levels for BitExtract
	logQ = append(logQ, 38, 38, 38)         // 3 levels for squaring
	logQ = append(logQ, 38)                 // extractExp
	logQ = append(logQ, 38)                 // Conv_{Real->Cplx}
	logQ = append(logQ, 38, 38, 38, 38, 38) // EvalCos
	logQ = append(logQ, 38)                 // Conv_{Cplx->Real}
	logQ = append(logQ, 38, 38, 38)         // CoeffsToSlots

	ciBase, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            logN,
		LogQ:            logQ,
		LogP:            []int{61, 61, 61, 61, 61, 61, 50},
		LogDefaultScale: 38,
		Xs:              ring.Ternary{H: 256},
		RingType:        ring.ConjugateInvariant,
	})
	if err != nil {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("ci params: %w", err)
	}
	stdLit := ciBase.ParametersLiteral()
	stdLit.RingType = ring.Standard
	params, err := ckks.NewParametersFromLiteral(stdLit)
	if err != nil {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("std params: %w", err)
	}

	// SlotsToCoeffs (homomorphic decoding).
	S2CParams := dft.MatrixLiteral{
		Type:     dft.HomomorphicDecode,
		LogSlots: params.LogMaxSlots(),
		LevelP:   params.MaxLevelP(),
		Levels:   []int{1},
	}
	S2CParams.LevelQ = len(S2CParams.Levels)

	// CoeffsToSlots (homomorphic encoding), real and imaginary parts split.
	C2SParams := dft.MatrixLiteral{
		Type:     dft.HomomorphicEncode,
		Format:   dft.SplitRealAndImag,
		LogSlots: params.LogMaxSlots(),
		LevelQ:   params.MaxLevel(),
		LevelP:   params.MaxLevelP(),
		Levels:   []int{1, 1, 1},
	}

	Mod1Params := mod1.ParametersLiteral{
		LevelQ:          params.MaxLevel() - C2SParams.Depth(true),
		LogScale:        38,
		Mod1Type:        mod1.CosDiscrete,
		Mod1Degree:      2 * ((1 << k) - 1),
		K:               1 << k,
		LogMessageRatio: k,
	}
	mod1P, err := mod1.NewParametersFromLiteral(params, Mod1Params)
	if err != nil {
		return ckks.Parameters{}, bootstrapping.Parameters{}, fmt.Errorf("mod1 params: %w", err)
	}
	F := (mod1P.ScalingFactor().Float64() / mod1P.MessageRatio()) / params.DefaultScale().Float64()
	S2CParams.Scaling = big.NewFloat(F)

	btpParams := bootstrapping.Parameters{
		ResidualParameters:      params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: S2CParams,
		Mod1ParametersLiteral:   Mod1Params,
		CoeffsToSlotsParameters: C2SParams,
		EphemeralSecretWeight:   32,
		CircuitOrder:            bootstrapping.DecodeThenModUp,
	}
	return params, btpParams, nil
}

// RefreshState refreshes the 128 bits of a CI state in groups of k via the CI->Std bridge,
// Algo1, BitExtract on both streams, recombination, and Std->CI return with a canonical scale.
func (c *Context) RefreshState(st aes.StateHE) (aes.StateHE, error) {
	eval := c.Eval

	// Flatten the 128 bits (order [byte][bit]).
	flat := make([]*rlwe.Ciphertext, 0, 128)
	for b := 0; b < 16; b++ {
		for i := 0; i < 8; i++ {
			flat = append(flat, st[b][i])
		}
	}

	// 1. CI -> Std bridge (z_j = A_j + i*B_j).
	stdFlat := make([]*rlwe.Ciphertext, len(flat))
	for i, ct := range flat {
		s, err := c.Sw.CIToStandard(ct)
		if err != nil {
			return aes.StateHE{}, fmt.Errorf("RefreshState CIToStandard %d: %w", i, err)
		}
		stdFlat[i] = s
	}

	// 2. Pack into 128/k integers m_A + i*m_B.
	groups, err := bitbatching.BitPackGroups(eval.Evaluator, stdFlat, c.K)
	if err != nil {
		return aes.StateHE{}, fmt.Errorf("RefreshState BitPackGroups: %w", err)
	}

	outFlat := make([]*rlwe.Ciphertext, len(flat))
	for g, ctPack := range groups {
		tPkt := time.Now()
		if d := ctPack.Level() - c.Target; d > 0 {
			eval.Evaluator.DropLevel(ctPack, d)
		}
		// 3. Algo1 (bootstrap) on the packed integers.
		ctReal, ctImag, err := algo1.Algo1(eval, c.Sw, ctPack, c.K)
		if err != nil {
			return aes.StateHE{}, fmt.Errorf("RefreshState Algo1 group %d: %w", g, err)
		}
		// 4. BitExtract each stream.
		bitsA, err := bitbatching.BitExtract(c.Params, eval.Evaluator, ctReal, c.K)
		if err != nil {
			return aes.StateHE{}, fmt.Errorf("RefreshState BitExtract A group %d: %w", g, err)
		}
		bitsB, err := bitbatching.BitExtract(c.Params, eval.Evaluator, ctImag, c.K)
		if err != nil {
			return aes.StateHE{}, fmt.Errorf("RefreshState BitExtract B group %d: %w", g, err)
		}

		// 5. Recombine + Std -> CI + canonical scale.
		for i := 0; i < c.K; i++ {
			z, err := utils.CombineReIm(eval.Evaluator, bitsA[i], bitsB[i])
			if err != nil {
				return aes.StateHE{}, fmt.Errorf("RefreshState combineReIm group %d bit %d: %w", g, i, err)
			}
			ci, err := c.Sw.StandardToCI(z)
			if err != nil {
				return aes.StateHE{}, fmt.Errorf("RefreshState StandardToCI group %d bit %d: %w", g, i, err)
			}
			// Option 2: FREE scale reset. Overwrite the scale metadata to the canonical 2^38
			// directly, instead of c.Sw.EvalCI.SetScale (which multiplies by the ratio then
			// Rescales -> costs a level). This is NOT value-preserving: the encrypted magnitude
			// and the tracked scale drift TOGETHER through the refresh Mul/Rescale chain, so
			// relabelling the scale transfers the accumulated deviation (old ci.Scale / c.Canon)
			// into the bit value. We accept it because (a) it puts all 128 refreshed bits on ONE
			// canonical scale, so the downstream MixColumns/AddRoundKey XOR-adds reconcile with no
			// scale-ratio loss, and (b) the round-final Cleaning (3x^2 - 2x^3) re-quantizes the
			// nudged bit back to 0/1. Only sound while old ci.Scale stays close to c.Canon: watch
			// the sc=2^... trace (now 3 decimals) to confirm the injected nudge is small.
			ci.Scale = c.Canon
			outFlat[g*c.K+i] = ci
		}
		fmt.Printf("   [refresh] packet %2d/%d done (%s)\n", g+1, len(groups), time.Since(tPkt).Round(time.Millisecond))
	}

	// Reassemble the state.
	var out aes.StateHE
	for b := 0; b < 16; b++ {
		for i := 0; i < 8; i++ {
			out[b][i] = outFlat[b*8+i]
		}
	}
	return out, nil
}

// AddRoundKeyCanon brings the round key to the state's (level, scale) with DropLevel + SetScale
// BEFORE the XOR. This is required: by AddRoundKey time the state is deep (low level) while the
// fresh round key sits at MaxLevel; the trick-mode xor aligns by squaring, which would square the
// key ~(MaxLevel - stateLevel) times in a row and blow up its encryption noise ((1+eps)^(2^n)).
// DropLevel descends the key without squaring, so no amplification. SetScale is applied only when
// the scale differs (it is an unconditional Mul+Rescale that would destroy the value at ratio
// 1.0), aiming one level above so the key lands exactly at the state level.
func (c *Context) AddRoundKeyCanon(st, rkHE aes.StateHE) (aes.StateHE, error) {
	L := st[0][0].Level()
	S := st[0][0].Scale
	var rkCanon aes.StateHE
	for b := 0; b < 16; b++ {
		for i := 0; i < 8; i++ {
			ct := rkHE[b][i].CopyNew()
			needScale := !ct.Scale.Equal(S)
			target := L
			if needScale {
				target = L + 1
			}
			if d := ct.Level() - target; d > 0 {
				c.Sw.EvalCI.DropLevel(ct, d)
			}
			if needScale {
				if err := c.Sw.EvalCI.SetScale(ct, S); err != nil {
					return aes.StateHE{}, fmt.Errorf("AddRoundKeyCanon SetScale byte %d bit %d: %w", b, i, err)
				}
			}
			rkCanon[b][i] = ct
		}
	}
	return c.AE.AddRoundKey(st, rkCanon)
}

// CleanState snaps every bit of the state to 0/1 (CI context).
func (c *Context) CleanState(st aes.StateHE) (aes.StateHE, error) {
	var out aes.StateHE
	for b := 0; b < 16; b++ {
		for i := 0; i < 8; i++ {
			cc, err := cleaning.Cleaning(c.Sw.CiP, c.Sw.EvalCI, st[b][i])
			if err != nil {
				return aes.StateHE{}, fmt.Errorf("CleanState byte %d bit %d: %w", b, i, err)
			}
			out[b][i] = cc
		}
	}
	return out, nil
}
