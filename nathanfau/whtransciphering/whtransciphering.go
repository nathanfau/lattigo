// Package whtransciphering runs AES under CKKS on the Walsh-Hadamard byte layout, refreshed by a
// conventional bootstrapping between rounds.
//
// One middle round is
//
//	four fused M_r' -> XOR tree -> AddRoundKey -> [Cleaning] -> Bootstrap
//
// and costs four levels: the M_r' are ct x pt, the XOR tree is two ct x ct products over the four
// T_r', AddRoundKey one more. On this encoding a XOR is a slot-wise multiplication, since
// WH(a) . WH(b) = WH(a xor b).
//
// Context.Round runs that sequence; each operation is exported on its own as well. The matrices
// come from nathanfau/whaes, the packing from blockpack.
package whtransciphering

import (
	"fmt"
	"runtime"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"

	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/blt"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockenc/charctx"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockpack"
	"github.com/tuneinsight/lattigo/v6/nathanfau/cleaning"
	"github.com/tuneinsight/lattigo/v6/nathanfau/params2"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/nathanfau/whaes"
)

// LastPattern names the fused matrix of the LAST round: no MixColumns there, so it is P_pi tensor
// A_1 alone. Fused takes it in place of an r' in 0..3.
const LastPattern = 4

// Config selects the variants a Context runs on. Clean has to be fixed here rather than flipped
// later: its depth decides how many primes the chain carries, so it is baked into the parameters.
type Config struct {
	Clean bool // (3x-x^3)/2 on the round output, two more levels
	Cache bool // keep the compiled M_r' resident: much faster, ~4x the RAM

	Chain params2.WHChain // the moduli chain and the bootstrapping knobs
}

type Context struct {
	Params    ckks.Parameters
	BtpParams bootstrapping.Parameters
	Ctx       *charctx.Context
	Ev        *blt.Evaluator
	Btp       *bootstrapping.Evaluator
	Cfg       Config

	Slots int
	SBox  [3][][]complex128 // A_1, A_2, A_3, shared by every fused matrix of every round

	RoundTop int // the level the refresh returns the state at, where a round starts

	// Levels the round keys have to be encrypted at.
	InitKeyLv int // rk0, just above the SlotsToCoeffs the first refresh opens on
	ARKKeyLv  int // middle rounds, after the XOR tree
	LastKeyLv int // last round, after the single fused matrix

	cache map[int]*blt.RawCompiled
}

// NewContext builds the whole machinery in its default configuration.
func NewContext(logN int) (*Context, error) {
	return NewContextWith(logN, Config{Chain: params2.DefaultWHChain()})
}

// NewContextWith is NewContext on the variants cfg selects. It generates the secret key, the
// relinearization key and the bootstrapping keys; the Galois keys of the fused matrices are added
// on demand by Fused, since they depend on the BSGS split of each one.
func NewContextWith(logN int, cfg Config) (*Context, error) {

	cfg.Chain.Clean = cfg.Clean // the chain and the circuit must agree on the depth
	params, btpParams, err := params2.WHTranscipheringParams(logN, cfg.Chain)
	if err != nil {
		return nil, fmt.Errorf("whtransciphering.NewContext: %w", err)
	}
	if params.MaxSlots() != whaes.Bytes*whaes.T {
		return nil, fmt.Errorf("whtransciphering.NewContext: %d slots, the layout needs %d",
			params.MaxSlots(), whaes.Bytes*whaes.T)
	}
	ctx, err := charctx.NewContext(params.ParametersLiteral())
	if err != nil {
		return nil, fmt.Errorf("whtransciphering.NewContext: charctx: %w", err)
	}

	// The round lives between SlotsToCoeffs and EvalMod, not at the top of the chain, so its
	// Galois keys are generated there. At the top they would be several times larger for levels
	// no fused matrix is ever evaluated at.
	roundTop := params2.WHRoundTop(cfg.Chain)
	ctx.GaloisLevelQ = &roundTop

	fmt.Println("BTS KeyGen ...")
	t0 := time.Now()
	btpKeys, _, err := btpParams.GenEvaluationKeys(ctx.SecretKey)
	if err != nil {
		return nil, fmt.Errorf("whtransciphering.NewContext: GenEvaluationKeys: %w", err)
	}
	btpEval, err := bootstrapping.NewEvaluator(btpParams, btpKeys)
	if err != nil {
		return nil, fmt.Errorf("whtransciphering.NewContext: bootstrapping evaluator: %w", err)
	}
	fmt.Printf("  bootstrapping keys in %s (%.0f MB)\n",
		time.Since(t0).Round(time.Millisecond), float64(btpKeys.BinarySize())/(1<<20))

	c := &Context{
		Params:    params,
		BtpParams: btpParams,
		Ctx:       ctx,
		Ev:        blt.NewEvaluatorWithCapacity(ctx, 1),
		Btp:       btpEval,
		Cfg:       cfg,
		Slots:     params.MaxSlots(),
		RoundTop:  roundTop,
		InitKeyLv: cfg.Chain.S2CDepth + 1,
		ARKKeyLv:  roundTop - 3,
		LastKeyLv: roundTop - 1,
		cache:     map[int]*blt.RawCompiled{},
	}

	for i, alpha := range whaes.Alphas {
		m, err := whaes.SBoxMatrix(alpha)
		if err != nil {
			return nil, fmt.Errorf("whtransciphering.NewContext: SBoxMatrix(%d): %w", alpha, err)
		}
		c.SBox[i] = m
	}
	return c, nil
}

// ---------------------------------------------------------------------------------------------
// encoding
// ---------------------------------------------------------------------------------------------

// Plaintext encodes a state at DefaultScale. The counter block is public, only the key is
// encrypted, so the initial AddRoundKey is a ct x pt product.
func (c *Context) Plaintext(s [whaes.Bytes]byte, level int) (*rlwe.Plaintext, error) {
	pt := ckks.NewPlaintext(c.Params, level)
	if err := c.Ctx.Encoder.Encode(blockpack.WHSlotVec(c.Slots, [][16]byte{s}), pt); err != nil {
		return nil, fmt.Errorf("whtransciphering.Plaintext: %w", err)
	}
	return pt, nil
}

// Encrypt is Plaintext under the secret key, for the round keys the client sends.
func (c *Context) Encrypt(s [whaes.Bytes]byte, level int) (*rlwe.Ciphertext, error) {
	pt, err := c.Plaintext(s, level)
	if err != nil {
		return nil, err
	}
	ct, err := c.Ctx.Encryptor.EncryptNew(pt)
	if err != nil {
		return nil, fmt.Errorf("whtransciphering.Encrypt: %w", err)
	}
	return ct, nil
}

// Decrypt reads the state back and returns the drift of the encoding, the worst | |value| - 1 |.
// No sign can have flipped while it stays below 1, so 1 - drift bounds the decision margin.
func (c *Context) Decrypt(ct *rlwe.Ciphertext) ([whaes.Bytes]byte, float64, error) {
	blocks, drift, err := blockpack.WHDecrypt(c.Params, c.Ctx.Encoder, c.Ctx.Decryptor, ct)
	if err != nil {
		return [whaes.Bytes]byte{}, 0, err
	}
	return blocks[0], drift, nil
}

// ---------------------------------------------------------------------------------------------
// the operations
// ---------------------------------------------------------------------------------------------

// compile encodes one fused matrix at the given input level. At LogN = 13 that is over a thousand
// diagonals of RingQP, several hundred megabytes, which is why they are built, used and released
// one at a time unless Cfg.Cache says otherwise.
func (c *Context) compile(pattern, level int) (*blt.RawCompiled, error) {
	if raw, ok := c.cache[pattern]; ok && raw.InputLevel == level {
		return raw, nil
	}

	var ks whaes.KronSum
	if pattern == LastPattern {
		ks = whaes.KronSum{
			Patterns: [][][]complex128{whaes.ShiftRowsMatrix()},
			Blocks:   [][][]complex128{c.SBox[0]},
		}
	} else {
		var err error
		if ks, err = whaes.RoundKronSum(pattern, c.SBox); err != nil {
			return nil, fmt.Errorf("whtransciphering.compile: %w", err)
		}
	}

	diagVals, diagonalOnly := ks.DiagonalValues(c.Slots)
	raw, err := blt.CompileDiagonalsWithOptions(
		lintrans.Diagonals[complex128](diagVals), diagonalOnly, nil,
		c.Slots, c.Slots, c.Ctx, level, blt.CompileOptions{})
	if err != nil {
		return nil, fmt.Errorf("whtransciphering.compile pattern %d: %w", pattern, err)
	}
	c.Ctx.EnsureGaloisKeys(raw.GaloisEls)

	if c.Cfg.Cache {
		c.cache[pattern] = raw
	}
	return raw, nil
}

// Fused applies M_r' for r' in 0..3, or P_pi tensor A_1 for LastPattern: the only ct x pt step of
// the round, one level. Without Cfg.Cache the compiled matrix is dropped inside the call, so the
// collector can take it back before the next one is built.
// It returns how long the compilation took, separately: encoding a thousand diagonals and applying
// them are different orders of magnitude, and only the first one Cfg.Cache can remove.
func (c *Context) Fused(st *rlwe.Ciphertext, pattern int) (*rlwe.Ciphertext, time.Duration, error) {
	t0 := time.Now()
	raw, err := c.compile(pattern, st.Level())
	if err != nil {
		return nil, 0, err
	}
	compile := time.Since(t0)

	out, err := c.Ev.ApplyRaw(st, raw)
	if err != nil {
		return nil, compile, fmt.Errorf("whtransciphering.Fused(%d): %w", pattern, err)
	}
	raw = nil
	if !c.Cfg.Cache {
		runtime.GC()
	}
	return out, compile, nil
}

// Xor is the whole point of the encoding: WH(a) . WH(b) = WH(a xor b) slot by slot, so a XOR is a
// single MulRelin and a Rescale. AddRoundKey is the same operation on a round key.
func (c *Context) Xor(a, b *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	x, y := a.CopyNew(), b.CopyNew()
	utils.AlignLevels(c.Ctx.Evaluator, x, y)
	out, err := c.Ctx.Evaluator.MulRelinNew(x, y)
	if err != nil {
		return nil, fmt.Errorf("whtransciphering.Xor: %w", err)
	}
	if err := c.Ctx.Evaluator.Rescale(out, out); err != nil {
		return nil, fmt.Errorf("whtransciphering.Xor rescale: %w", err)
	}
	return out, nil
}

// AddRoundKey is Xor on the encrypted round key.
func (c *Context) AddRoundKey(st, rk *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return c.Xor(st, rk)
}

// Clean pushes the slots back onto +-1 with (3x - x^3)/2, and costs the two extra primes of
// Cfg.Clean.
func (c *Context) Clean(st *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	return cleaning.SignCleaning(c.Params, c.Ctx.Evaluator, st)
}

// Bootstrap refreshes the state, assembled step by step rather than through btp.Bootstrap
func (c *Context) Bootstrap(st *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	ct, err := c.Btp.SlotsToCoeffs(st, nil)
	if err != nil {
		return nil, fmt.Errorf("whtransciphering.Bootstrap SlotsToCoeffs: %w", err)
	}

	if ct, _, err = c.Btp.ScaleDown(ct); err != nil {
		return nil, fmt.Errorf("whtransciphering.Bootstrap ScaleDown: %w", err)
	}

	if ct, err = c.Btp.ModUp(ct); err != nil {
		return nil, fmt.Errorf("whtransciphering.Bootstrap ModUp: %w", err)
	}

	ctReal, ctImag, err := c.Btp.CoeffsToSlots(ct)
	if err != nil {
		return nil, fmt.Errorf("whtransciphering.Bootstrap CoeffsToSlots: %w", err)
	}
	if ctImag == nil {
		return nil, fmt.Errorf("whtransciphering.Bootstrap: CoeffsToSlots returned one half, the layout needs both")
	}

	if ctReal, err = c.Btp.EvalMod(ctReal); err != nil {
		return nil, fmt.Errorf("whtransciphering.Bootstrap EvalMod real: %w", err)
	}
	if ctImag, err = c.Btp.EvalMod(ctImag); err != nil {
		return nil, fmt.Errorf("whtransciphering.Bootstrap EvalMod imag: %w", err)
	}

	out, err := utils.CombineReIm(c.Ctx.Evaluator, ctReal, ctImag)
	if err != nil {
		return nil, fmt.Errorf("whtransciphering.Bootstrap CombineReIm: %w", err)
	}
	out.Scale = c.Params.DefaultScale()
	return out, nil
}

// ---------------------------------------------------------------------------------------------
// the rounds
// ---------------------------------------------------------------------------------------------

// RoundStep is called after each operation: its name, its duration, the state it produced. That is
// where a caller hooks its timings and its oracle checks; nil skips it.
//
// A nil ciphertext means the call carries a duration only, with no state to check: it is how the
// breakdown of an operation is reported next to the operation itself, so a caller can record both
// without decrypting twice.
type RoundStep func(name string, dur time.Duration, ct *rlwe.Ciphertext)

func noStep(string, time.Duration, *rlwe.Ciphertext) {}

// FirstRound is AES round 0, the initial AddRoundKey. The block is public, so it goes in as a
// plaintext and the product costs one level. rk0 comes in at InitKeyLv so that the product lands
// exactly where the refresh opens its SlotsToCoeffs; any level spent above would be thrown away
// by the ScaleDown that follows.
func (c *Context) FirstRound(block [whaes.Bytes]byte, rk0 *rlwe.Ciphertext, after RoundStep) (*rlwe.Ciphertext, error) {
	if after == nil {
		after = noStep
	}

	pt, err := c.Plaintext(block, rk0.Level())
	if err != nil {
		return nil, fmt.Errorf("FirstRound: %w", err)
	}

	t0 := time.Now()
	st, err := c.Ctx.Evaluator.MulNew(rk0, pt)
	if err != nil {
		return nil, fmt.Errorf("FirstRound: %w", err)
	}
	if err := c.Ctx.Evaluator.Rescale(st, st); err != nil {
		return nil, fmt.Errorf("FirstRound rescale: %w", err)
	}
	after("AddRoundKey", time.Since(t0), st)

	t0 = time.Now()
	if st, err = c.Bootstrap(st); err != nil {
		return nil, fmt.Errorf("FirstRound Bootstrap: %w", err)
	}
	after("Bootstrap", time.Since(t0), st)
	return st, nil
}

// Round advances the state by one AES middle round. rk is the round key at ARKKeyLv. The state
// enters and leaves at the top of the chain, the bootstrap being the last step.
//
// after sees the four T_r' under the names M_0..M_3, the two halves of the XOR tree, then
// "XOR tree", "AddRoundKey", "Cleaning" and "Bootstrap". Without Cfg.Cache, the M_r' durations
// include the compilation of the matrix, which is what the round really pays.
func (c *Context) Round(st, rk *rlwe.Ciphertext, after RoundStep) (*rlwe.Ciphertext, error) {
	if after == nil {
		after = noStep
	}

	var Ts [4]*rlwe.Ciphertext
	for rPrime := 0; rPrime < 4; rPrime++ {
		t0 := time.Now()
		out, compile, err := c.Fused(st, rPrime)
		if err != nil {
			return nil, fmt.Errorf("Round: %w", err)
		}
		Ts[rPrime] = out
		total := time.Since(t0)
		after(fmt.Sprintf("M_%d compile", rPrime), compile, nil)
		after(fmt.Sprintf("M_%d apply", rPrime), total-compile, nil)
		after(fmt.Sprintf("M_%d", rPrime), total, out)
	}

	t0 := time.Now()
	left, err := c.Xor(Ts[0], Ts[1])
	if err != nil {
		return nil, fmt.Errorf("Round T0^T1: %w", err)
	}
	after("T0^T1", time.Since(t0), left)

	t0 = time.Now()
	right, err := c.Xor(Ts[2], Ts[3])
	if err != nil {
		return nil, fmt.Errorf("Round T2^T3: %w", err)
	}
	after("T2^T3", time.Since(t0), right)

	t0 = time.Now()
	if st, err = c.Xor(left, right); err != nil {
		return nil, fmt.Errorf("Round XOR tree: %w", err)
	}
	after("XOR tree", time.Since(t0), st)

	return c.arkThenRefresh(st, rk, after)
}

// LastRound is AES round 10: no MixColumns
func (c *Context) LastRound(st, rk *rlwe.Ciphertext, after RoundStep) (*rlwe.Ciphertext, error) {
	if after == nil {
		after = noStep
	}

	t0 := time.Now()
	st, compile, err := c.Fused(st, LastPattern)
	if err != nil {
		return nil, fmt.Errorf("LastRound: %w", err)
	}
	total := time.Since(t0)
	after("SubBytes+SR compile", compile, nil)
	after("SubBytes+SR apply", total-compile, nil)
	after("SubBytes+SR", total, st)

	t0 = time.Now()
	if st, err = c.AddRoundKey(st, rk); err != nil {
		return nil, fmt.Errorf("LastRound AddRoundKey: %w", err)
	}
	after("AddRoundKey", time.Since(t0), st)

	t0 = time.Now()
	if st, err = c.Clean(st); err != nil {
		return nil, fmt.Errorf("LastRound Cleaning: %w", err)
	}
	after("Cleaning", time.Since(t0), st)
	return st, nil
}

// arkThenRefresh closes a middle round: AddRoundKey, the optional cleaning, then the bootstrap.
func (c *Context) arkThenRefresh(st, rk *rlwe.Ciphertext, after RoundStep) (*rlwe.Ciphertext, error) {
	t0 := time.Now()
	st, err := c.AddRoundKey(st, rk)
	if err != nil {
		return nil, fmt.Errorf("AddRoundKey: %w", err)
	}
	after("AddRoundKey", time.Since(t0), st)

	if c.Cfg.Clean {
		t0 = time.Now()
		if st, err = c.Clean(st); err != nil {
			return nil, fmt.Errorf("Cleaning: %w", err)
		}
		after("Cleaning", time.Since(t0), st)
	}

	t0 = time.Now()
	if st, err = c.Bootstrap(st); err != nil {
		return nil, fmt.Errorf("Bootstrap: %w", err)
	}
	after("Bootstrap", time.Since(t0), st)
	return st, nil
}
