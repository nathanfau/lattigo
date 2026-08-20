package whtransciphering

//	go test ./nathanfau/whtransciphering/ -run '^TestWHAES$' -v -timeout 0

import (
	"flag"
	"fmt"
	"math/rand"
	rtdebug "runtime/debug" // nathanfau/debug already owns the name
	"strings"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"

	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockpack"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/nathanfau/params2"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/nathanfau/whaes"
)

var (
	roundsFlag = flag.Int("rounds", 9, "number of AES middle rounds; 9 is the full AES-128")
	seedFlag   = flag.Int64("seed", 1, "seed of the block draw")
)

// reclaim hands a fused matrix back to the OS the moment it has been applied.
//
// One M_r' is 1277 diagonals of (round level + 1) + (levelP + 1) polynomials, ~1.09 GB at logN 13,
// and a round builds four of them in a row. Fused drops each one as soon as it is applied, but the
// collector only acts when its own heap target is reached, so up to four dead matrices pile up:
// measured, the peak goes to 2.9 GB and the fourth matrix takes 32 s to encode against 5 s for the
// first, the run being by then busy collecting rather than computing. Returning each one keeps the
// round to the single matrix it actually needs at a time. Same argument, more so, for the
// SubBytes+SR of the last round, which carries 2036 diagonals.
func reclaim(step string) {
	if strings.HasPrefix(step, "M_") || step == "SubBytes+SR" {
		rtdebug.FreeOSMemory()
	}
}

// config reads the flags
func config() Config {
	return Config{Clean: true, Chain: params2.DefaultWHChain()}
}

// ---------------------------------------------------------------------------------------------
// the cleartext oracle
// ---------------------------------------------------------------------------------------------

// tState is the byte state T_r' carries. M_r' sends output block 4c+r to A_{alpha_{r,r'}} applied
// to input block pi^{-1}(4c+r'), and A_alpha sends WH(m) to WH(alpha . S(m)), so T_r' is a genuine
// WH encoding of a state -- one that can be decoded and compared byte by byte.
func tState(s [whaes.Bytes]byte, rPrime int) [whaes.Bytes]byte {
	sub := s
	aes.SubBytes(sub[:])
	aes.ShiftRows(sub[:])

	var out [whaes.Bytes]byte
	for c := 0; c < 4; c++ {
		for r := 0; r < 4; r++ {
			out[4*c+r] = aes.GFMul(whaes.MixAlpha[r][rPrime], sub[4*c+rPrime])
		}
	}
	return out
}

// xorState is what the slot-wise product of two WH states decodes to.
func xorState(a, b [whaes.Bytes]byte) (out [whaes.Bytes]byte) {
	for i := range out {
		out[i] = a[i] ^ b[i]
	}
	return out
}

// diffBytes counts how many bytes differ and returns the first differing index.
func diffBytes(got, want [whaes.Bytes]byte) (wrong, first int) {
	first = -1
	for b := 0; b < whaes.Bytes; b++ {
		if got[b] != want[b] {
			if first < 0 {
				first = b
			}
			wrong++
		}
	}
	return wrong, first
}

// stepWant is the state each named step of a round should hold, given the state at round entry.
// The round advances the oracle only once, at "XOR tree"; everything before is an intermediate.
func stepWant(name string, entry, mixed, keyed [whaes.Bytes]byte) [whaes.Bytes]byte {
	switch name {
	case "M_0", "M_1", "M_2", "M_3":
		return tState(entry, int(name[2]-'0'))
	case "T0^T1":
		return xorState(tState(entry, 0), tState(entry, 1))
	case "T2^T3":
		return xorState(tState(entry, 2), tState(entry, 3))
	case "XOR tree":
		return mixed
	default: // AddRoundKey, Cleaning, Bootstrap
		return keyed
	}
}

// ---------------------------------------------------------------------------------------------
// reporting
// ---------------------------------------------------------------------------------------------

// report prints the slots, the chain and the precision of one step, then the byte oracle. It is
// the counterpart of transciphering's own report, on a single ciphertext.
func report(c *Context, round int, step string, ct *rlwe.Ciphertext, want [whaes.Bytes]byte) {
	label := fmt.Sprintf("T%d %-11s", round, step)
	debug.DbgSlotStd(label+" ct    =", ct)
	debug.DbgChain(label+" chain :", c.Ctx.Evaluator, ct)
	debug.PrecPoolStd(label+" prec  :",
		debug.PrecCt{Name: step, Want: blockpack.WHSlotVec(c.Slots, [][16]byte{want}), Ct: ct})

	got, drift, err := c.Decrypt(ct)
	if err != nil {
		fmt.Printf("  [T%d] %-11s ORACLE: decode error: %v\n", round, step, err)
		return
	}
	if got == want {
		fmt.Printf("  [T%d] %-11s ORACLE: TRUE  %x  (drift %.4f, sign margin >= %.4f)\n",
			round, step, got, drift, 1-drift)
		return
	}
	wrong, first := diffBytes(got, want)
	fmt.Printf("  [T%d] %-11s ORACLE: FALSE %d/16 bytes, first byte %d (drift %.4f)\n            got =%x\n            want=%x\n",
		round, step, wrong, first, drift, got, want)
}

// printChain cuts the one moduli chain into the stages that consume it. The circuit runs top-down
// and loops back to q0, so read from the bottom: q0, SlotsToCoeffs, the AES round, EvalMod, then
// CoeffsToSlots at the top.
func printChain(c *Context, btp bootstrapping.Parameters) {
	sum := func(v []int) (s int) {
		for _, x := range v {
			s += x
		}
		return s
	}
	full := btp.BootstrappingParameters.LogQi()
	nS2C := sum(btp.SlotsToCoeffsParameters.Levels)
	nRound := params2.RoundDepth(c.Cfg.Chain)
	nMod1 := btp.Mod1ParametersLiteral.Depth()
	nC2S := sum(btp.CoeffsToSlotsParameters.Levels)

	if 1+nS2C+nRound+nMod1+nC2S != len(full) {
		fmt.Printf("  chain            : %d  (round %d, S2C %d, EvalMod %d, C2S %d)\n",
			full, nRound, nS2C, nMod1, nC2S)
		return
	}

	seg := func(name string, from, n int) {
		if n == 0 {
			return
		}
		fmt.Printf("  %-16s levels %2d..%-2d  %d\n", name, from, from+n-1, full[from:from+n])
	}
	seg("q0", 0, 1)
	seg("SlotsToCoeffs", 1, nS2C)
	seg("AES round", 1+nS2C, nRound)
	seg("EvalMod", 1+nS2C+nRound, nMod1)
	seg("CoeffsToSlots", 1+nS2C+nRound+nMod1, nC2S)

	fmt.Printf("  %-16s %d levels = C2S %d + EvalMod %d + S2C %d\n",
		"refresh depth", nS2C+nMod1+nC2S, nC2S, nMod1, nS2C)
	fmt.Printf("  %-16s K=%d, degree %d, %d double angles, message ratio 2^%d, H(%d dense; %d ephemeral)\n",
		"EvalMod",
		btp.Mod1ParametersLiteral.K,
		btp.Mod1ParametersLiteral.Mod1Degree,
		btp.Mod1ParametersLiteral.DoubleAngle,
		btp.Mod1ParametersLiteral.LogMessageRatio,
		btp.BootstrappingParameters.XsHammingWeight(),
		btp.EphemeralSecretWeight)
}

// ---------------------------------------------------------------------------------------------
// the test
// ---------------------------------------------------------------------------------------------

// TestWHAES runs the whole cipher on one block: FirstRound -> n * Round -> LastRound, with the
// oracle checked after EVERY operation.
func TestWHAES(t *testing.T) {

	if testing.Short() {
		t.Skip("full AES over WH blocks: a bootstrapping per round, several minutes")
	}

	nMiddle := *roundsFlag
	if nMiddle < 0 || nMiddle > 9 {
		t.Fatalf("-rounds %d, want 0..9", nMiddle)
	}
	last := nMiddle + 1

	c, err := NewContextWith(13, config())
	if err != nil {
		t.Fatalf("NewContext: %v", err)
	}
	debug.EncStd, debug.DecStd, debug.ParamsStd = c.Ctx.Encoder, c.Ctx.Decryptor, c.Params

	fmt.Println("================ parameters ================")
	debug.DbgParams("WH", c.BtpParams.BootstrappingParameters)
	printChain(c, c.BtpParams)
	fmt.Printf("  %-16s %d levels, %d..%d (BLT 1 + XOR tree 2 + ARK 1 + cleaning 2)\n", "round budget",
		params2.RoundDepth(c.Cfg.Chain), c.RoundTop, c.RoundTop-params2.RoundDepth(c.Cfg.Chain)+1)
	fmt.Printf("  %-16s rk0 %d, middle %d, last %d\n", "round keys", c.InitKeyLv, c.ARKKeyLv, c.LastKeyLv)

	// ---- AES material ----

	key := [16]byte{0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c}
	rk := aes.KeyExpansion(key[:])

	rng := rand.New(rand.NewSource(*seedFlag))
	var block [whaes.Bytes]byte
	for i := range block {
		block[i] = byte(rng.Intn(256))
	}

	fmt.Printf("\n AES-128: %d middle rounds + last, one 16-byte block over %d WH slots, seed = %d \n",
		nMiddle, c.Slots, *seedFlag)

	// The round keys come from the client, encrypted at the level where they meet the state.
	rkCT := make([]*rlwe.Ciphertext, len(rk))
	encRK := func(r, level int) *rlwe.Ciphertext {
		var s [whaes.Bytes]byte
		copy(s[:], rk[r])
		ct, err := c.Encrypt(s, level)
		if err != nil {
			t.Fatalf("encrypt round key %d: %v", r, err)
		}
		return ct
	}
	rkCT[0] = encRK(0, c.InitKeyLv)
	for r := 1; r <= nMiddle; r++ {
		rkCT[r] = encRK(r, c.ARKKeyLv)
	}
	rkCT[last] = encRK(last, c.LastKeyLv)

	// ---- the run ----

	state := block
	times := map[string][]time.Duration{}
	round := 0
	entry, mixed, keyed := state, state, state

	after := func(name string, dur time.Duration, ct *rlwe.Ciphertext) {
		times[name] = append(times[name], dur)
		if ct == nil { // a timing breakdown, not a state: nothing to check
			return
		}
		reclaim(name)
		report(c, round, name, ct, stepWant(name, entry, mixed, keyed))
	}

	tGlobal := time.Now()

	fmt.Println("\n================ Round 0 (initial ARK, ct x pt) ================")
	aes.AddRoundKey(state[:], rk[0])
	entry, mixed, keyed = state, state, state
	st, err := c.FirstRound(block, rkCT[0], after)
	if err != nil {
		t.Fatalf("FirstRound: %v", err)
	}

	for r := 1; r <= nMiddle; r++ {
		fmt.Printf("\n================ Round %d/%d ================\n", r, last)
		round = r

		// The oracle of the whole round, computed before it runs: stepWant picks from it.
		entry = state
		mixed = state
		aes.SubBytes(mixed[:])
		aes.ShiftRows(mixed[:])
		aes.MixColumns(mixed[:])
		keyed = mixed
		aes.AddRoundKey(keyed[:], rk[r])

		tRound := time.Now()
		if st, err = c.Round(st, rkCT[r], after); err != nil {
			t.Fatalf("Round %d: %v", r, err)
		}
		times["Round"] = append(times["Round"], time.Since(tRound))
		fmt.Printf("  [T%d] round time (incl. oracle traces): %s\n", r, time.Since(tRound).Round(time.Millisecond))

		state = keyed
	}

	fmt.Printf("\n================ Round %d/%d (last, no MixColumns) ================\n", last, last)
	round = last
	entry = state
	mixed = state
	aes.SubBytes(mixed[:])
	aes.ShiftRows(mixed[:])
	keyed = mixed
	aes.AddRoundKey(keyed[:], rk[last])

	tRound := time.Now()
	// "SubBytes+SR" stops at mixed, AddRoundKey lands on keyed: stepWant's default covers the
	// second, the first needs the intermediate.
	lastAfter := func(name string, dur time.Duration, ct *rlwe.Ciphertext) {
		times[name] = append(times[name], dur)
		if ct == nil {
			return
		}
		want := keyed
		if name == "SubBytes+SR" {
			want = mixed
		}
		reclaim(name)
		report(c, round, name, ct, want)
	}
	if st, err = c.LastRound(st, rkCT[last], lastAfter); err != nil {
		t.Fatalf("LastRound: %v", err)
	}
	times["LastRound"] = append(times["LastRound"], time.Since(tRound))
	state = keyed

	// ---- result ----

	fmt.Printf("\nTOTAL (incl. oracle traces): %s\n", time.Since(tGlobal).Round(time.Millisecond))
	fmt.Println("================ timing stats ================")
	ops := []string{"M_0", "M_1", "M_2", "M_3", "T0^T1", "T2^T3", "XOR tree", "AddRoundKey", "Cleaning", "Bootstrap"}
	var opMeanSum time.Duration
	for _, name := range ops {
		if ds := times[name]; len(ds) > 0 {
			var s time.Duration
			for _, d := range ds {
				s += d
			}
			opMeanSum += s / time.Duration(len(ds))
		}
	}
	for _, name := range ops {
		utils.TimeStats(name, times[name], opMeanSum)
	}
	utils.TimeStats("Round", times["Round"])
	utils.TimeStats("LastRound", times["LastRound"])

	// The two halves of a fused matrix: encoding its diagonals, then applying them. Only the first
	// is what Cfg.Cache removes, so the split says whether caching is worth its memory. Shares are
	// against the same round total as above, so they read against the M_i lines.
	fmt.Println("---------------- of which, per fused matrix ----------------")
	for _, name := range []string{
		"M_0 compile", "M_0 apply", "M_1 compile", "M_1 apply",
		"M_2 compile", "M_2 apply", "M_3 compile", "M_3 apply",
		"SubBytes+SR compile", "SubBytes+SR apply",
	} {
		utils.TimeStats(name, times[name], opMeanSum)
	}

	got, drift, err := c.Decrypt(st)
	if err != nil {
		t.Fatalf("final decode: %v", err)
	}
	fmt.Printf("\n  block  : %x\n  got    : %x\n  want   : %x\n  drift: %.4f (sign margin >= %.4f)\n",
		block, got, state, drift, 1-drift)

	if got != state {
		wrong, first := diffBytes(got, state)
		t.Fatalf("WH AES: %d/16 bytes wrong, first byte %d", wrong, first)
	}
	if nMiddle == 9 {
		var ref [whaes.Bytes]byte
		copy(ref[:], aes.EncryptBlock(block[:], key[:]))
		if got != ref {
			t.Fatalf("matches the traced oracle but not FIPS-197: got %x, want %x", got, ref)
		}
		fmt.Printf("\n=== OK: full AES-128, the block matches FIPS-197 on all %d bytes ===\n", whaes.Bytes)
		return
	}
	fmt.Printf("\n=== OK: %d middle round(s) + last, all %d bytes conform to the AES oracle ===\n", nMiddle, whaes.Bytes)
}
