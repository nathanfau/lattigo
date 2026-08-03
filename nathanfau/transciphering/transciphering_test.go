package transciphering

//	go test ./nathanfau/transciphering/ -run '^TestTransciphering$' -v -subbytes 2 -rounds 1 -timeout 0
//	go test ./nathanfau/transciphering/ -run '^TestAES$' -v -subbytes 2 -timeout 0

import (
	"flag"
	"fmt"
	"math/rand"

	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockpack"
	"github.com/tuneinsight/lattigo/v6/nathanfau/cleaning"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
)

var (
	nRoundsFlag = flag.Int("rounds", 1, "number of AES middle rounds to run")
	sbVersion   = flag.Int("subbytes", 2, "SubBytes version (1..3), by decreasing cost: 247, 98, 69 relin per byte")
	xorFlag     = flag.String("xor", "nosq", `XOR circuit used throughout: "nosq" for x+y-2xy, "sq" for (x-y)^2`)
	cleanFlag   = flag.String("clean", "cleaning", `cleaning polynomial: "cleaning" (2 levels), "smoother" or "verysmoother" (3 levels, one prime more)`)
	placeFlag   = flag.String("place", "after", `where the round cleaning lands: "after" (AddRoundKey then Cleaning), "both" (clean state and key, then XOR) or "one" (clean the state only, then XOR)`)
	seedFlag    = flag.Int64("seed", 0, "seed of the block draw; 0 draws one from the clock")
)

// blockSeed fixes WHICH blocks the batch carries, so two configs can be compared on the same
// input. It does NOT fix the key or the encryption noise, which lattigo draws from crypto/rand:
// runs on one seed are far closer than on two, not identical.
func blockSeed() int64 {
	if *seedFlag != 0 {
		return *seedFlag
	}
	return time.Now().UnixNano()
}

// config parses the flags before any keygen, so a typo fails in a second instead of a minute.
func config(t *testing.T) Config {
	t.Helper()
	xk, err := aes.ParseXorKind(*xorFlag)
	if err != nil {
		t.Fatalf("-xor: %v", err)
	}
	ck, err := cleaning.ParseKind(*cleanFlag)
	if err != nil {
		t.Fatalf("-clean: %v", err)
	}
	pk, err := cleaning.ParsePlacement(*placeFlag)
	if err != nil {
		t.Fatalf("-place: %v", err)
	}
	return Config{Xor: xk, Clean: ck, Place: pk}
}

// tail names the steps that close a round, which every placement but "after" merges into one
func tail(cfg Config) []string {
	if cfg.Place != cleaning.CleanAfter {
		return []string{"XorClean"}
	}
	return []string{"AddRoundKey", "Cleaning"}
}

// TestTransciphering runs Context.Round on a full batch of DISTINCT random blocks, one per slot,
// with a slot / chain / precision trace and a row-major AES-oracle comparison after EVERY
// operation. The packed layout is invariant, so rounds > 1 chain with no return trip.
func TestTransciphering(t *testing.T) {
	const logN, k = 11, 4

	n := *nRoundsFlag
	if n < 1 {
		n = 1
	}
	cfg := config(t)
	fmt.Printf(" Transciphering: rounds=%d, SubBytes=V%d, XOR=%s, clean=%s, place=%s \n",
		n, *sbVersion, cfg.Xor, cfg.Clean, cfg.Place)

	ctx, err := NewContextWith(logN, k, cfg)
	if err != nil {
		t.Fatalf("NewContext: %v", err)
	}
	debug.DbgParams("TranscipheringParams", ctx.Params)
	ciP := ctx.Sw.CiP

	// FIPS-197 key (shared key stream). Each slot carries a DISTINCT random block (real batching);
	// states[bi] tracks block bi's cleartext state, and entry to round 1 is plaintext XOR rk0.
	key := [16]byte{0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c}
	rk := aes.KeyExpansion(key[:])

	seed := blockSeed()
	states := blockpack.RandomBlocks(ciP, rand.New(rand.NewSource(seed)))
	fmt.Printf("random seed = %d  (%d blocks)\n", seed, len(states))
	for bi := range states {
		aes.AddRoundKey(states[bi][:], rk[0])
	}

	st, err := blockpack.Encrypt(ciP, ctx.EcdCI, ctx.EncCI, states, SubBytesLevel)
	if err != nil {
		t.Fatalf("blockpack.Encrypt state: %v", err)
	}

	rkHE := make([]blockpack.Packed, n+1)
	for r := 1; r <= n; r++ {
		rkHE[r] = encRK(t, ctx, rk[r], len(states), ctx.ARKKeyLv)
	}

	fmt.Println("================ input (entry to round 1) ================")
	report(ctx, st, states, 0, "input")

	times := map[string][]time.Duration{}
	tGlobal := time.Now()

	for r := 1; r <= n; r++ {
		fmt.Printf("\n================ Round %d/%d ================\n", r, n)
		tRound := time.Now()

		// after times the operation, advances the cleartext oracle by the same one and prints the
		// trace. Refresh advances it by ShiftRows: the bootstrap keeps the bit values, only the
		// ShiftRows at the Algo1 pause moves them.
		after := func(name string, dur time.Duration, out blockpack.Packed) {
			st = out
			times[name] = append(times[name], dur)
			for bi := range states {
				switch name {
				case "SubBytes":
					aes.SubBytes(states[bi][:])
				case "Refresh":
					states[bi] = aes.ShiftRowsRM(states[bi])
				case "MixColumns":
					states[bi] = aes.MixColumnsRM(states[bi])
				case "AddRoundKey", "XorClean":
					aes.AddRoundKey(states[bi][:], rk[r])
				}
			}
			report(ctx, st, states, r, name)
		}

		if st, err = ctx.Round(st, rkHE[r], *sbVersion, after); err != nil {
			t.Fatalf("Round T%d: %v", r, err)
		}

		dRound := time.Since(tRound)
		times["Round"] = append(times["Round"], dRound)
		fmt.Printf("  [T%d] round time (incl. oracle traces): %s\n", r, dRound.Round(time.Millisecond))
	}

	fmt.Printf("\nTOTAL (%d rounds, incl. oracle traces): %s\n", n, time.Since(tGlobal).Round(time.Millisecond))
	fmt.Println("================ timing stats (per round) ================")
	// Share reference = sum of the op means, so the ops add up to ~100%; "Round" is the wall time
	// (ops + oracle traces) and gets no share.
	ops := append([]string{"SubBytes", "Refresh", "MixColumns"}, tail(cfg)...)
	var opMeanSum time.Duration
	for _, name := range ops {
		if ds := times[name]; len(ds) > 0 {
			var sum time.Duration
			for _, d := range ds {
				sum += d
			}
			opMeanSum += sum / time.Duration(len(ds))
		}
	}
	for _, name := range ops {
		utils.TimeStats(name, times[name], opMeanSum)
	}
	utils.TimeStats("Round", times["Round"])

	if err := blockpack.Check(ciP, ctx.EcdCI, ctx.DecCI, st, states); err != nil {
		t.Errorf("AES (%d rounds): %v", n, err)
	} else {
		fmt.Printf("\n=== OK: %d middle round(s), all %d blocks conform to the row-major AES oracle ===\n", n, len(states))
	}
}

// TestAES runs the whole cipher on a full batch of DISTINCT random blocks:
// FirstRound -> 9 * Round -> LastRoundV1, with the oracle checked after EVERY operation.
func TestAES(t *testing.T) {
	const logN, k = 12, 4

	if testing.Short() {
		t.Skip("full AES: 10 rounds, several minutes")
	}

	cfg := config(t)

	ctx, err := NewContextWith(logN, k, cfg)
	if err != nil {
		t.Fatalf("NewContext: %v", err)
	}
	debug.DbgParams("TranscipheringParams", ctx.Params)
	ciP := ctx.Sw.CiP

	key := [16]byte{0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c}
	rk := aes.KeyExpansion(key[:])
	nMiddle := len(rk) - 2 // round 0 is the initial ARK, the last one has no MixColumns

	seed := blockSeed()
	blocks := blockpack.RandomBlocks(ciP, rand.New(rand.NewSource(seed)))
	fmt.Printf(" AES-128: %d middle rounds, SubBytes=V%d, XOR=%s, clean=%s, place=%s, %d blocks, random seed = %d \n",
		nMiddle, *sbVersion, cfg.Xor, cfg.Clean, cfg.Place, len(blocks), seed)

	// The blocks stay in the CLEAR: that is what the client sends and what FirstRound takes.
	states := make([][16]byte, len(blocks))
	copy(states, blocks)

	// Three key levels: rk0 lands before SubBytes, the middle keys after MixColumns, the last one
	// straight out of the refresh since no MixColumns eats into it.
	rk0 := encRK(t, ctx, rk[0], len(blocks), InitLevel)
	rkHE := make([]blockpack.Packed, len(rk))
	for r := 1; r <= nMiddle; r++ {
		rkHE[r] = encRK(t, ctx, rk[r], len(blocks), ctx.ARKKeyLv)
	}
	rkHE[len(rk)-1] = encRK(t, ctx, rk[len(rk)-1], len(blocks), ctx.LastKeyLv)

	times := map[string][]time.Duration{}
	tGlobal := time.Now()

	fmt.Println("\n================ Round 0 (initial ARK) ================")
	t0 := time.Now()
	st, err := ctx.FirstRound(blocks, rk0)
	if err != nil {
		t.Fatalf("FirstRound: %v", err)
	}
	times["FirstRound"] = append(times["FirstRound"], time.Since(t0))
	for bi := range states {
		aes.AddRoundKey(states[bi][:], rk[0])
	}
	if l := st[0][0].Level(); l != SubBytesLevel {
		t.Errorf("FirstRound left the state at level %d, want SubBytesLevel %d", l, SubBytesLevel)
	}
	report(ctx, st, states, 0, "FirstRound")

	// round and rkNow name the round in flight, so the same hook serves the middle rounds and the
	// last one. Refresh advances the oracle by ShiftRows: the bootstrap keeps the bit values, only
	// the ShiftRows at the Algo1 pause moves them.
	round, rkNow := 0, rk[0]
	after := func(name string, dur time.Duration, out blockpack.Packed) {
		st = out
		times[name] = append(times[name], dur)
		for bi := range states {
			switch name {
			case "SubBytes":
				aes.SubBytes(states[bi][:])
			case "Refresh":
				states[bi] = aes.ShiftRowsRM(states[bi])
			case "MixColumns":
				states[bi] = aes.MixColumnsRM(states[bi])
			case "AddRoundKey", "XorClean":
				aes.AddRoundKey(states[bi][:], rkNow)
			}
		}
		report(ctx, st, states, round, name)
	}

	for r := 1; r <= nMiddle; r++ {
		fmt.Printf("\n================ Round %d/%d ================\n", r, len(rk)-1)
		round, rkNow = r, rk[r]
		tRound := time.Now()
		if st, err = ctx.Round(st, rkHE[r], *sbVersion, after); err != nil {
			t.Fatalf("Round T%d: %v", r, err)
		}
		times["Round"] = append(times["Round"], time.Since(tRound))
	}

	last := len(rk) - 1
	fmt.Printf("\n================ Round %d/%d (last, no MixColumns) ================\n", last, last)
	round, rkNow = last, rk[last]
	tRound := time.Now()
	if st, err = ctx.LastRoundV1(st, rkHE[last], *sbVersion, after); err != nil {
		t.Fatalf("LastRoundV1: %v", err)
	}
	times["LastRound"] = append(times["LastRound"], time.Since(tRound))

	fmt.Printf("\nTOTAL AES-128 (incl. oracle traces): %s\n", time.Since(tGlobal).Round(time.Millisecond))
	fmt.Println("================ timing stats ================")
	// Share reference = sum of the op means, so the ops add up to ~100%; the round totals are wall
	// time (ops + oracle traces) and get no share.
	ops := append([]string{"SubBytes", "Refresh", "MixColumns"}, tail(cfg)...)
	var opMeanSum time.Duration
	for _, name := range ops {
		if ds := times[name]; len(ds) > 0 {
			var sum time.Duration
			for _, d := range ds {
				sum += d
			}
			opMeanSum += sum / time.Duration(len(ds))
		}
	}
	for _, name := range ops {
		utils.TimeStats(name, times[name], opMeanSum)
	}
	utils.TimeStats("FirstRound", times["FirstRound"])
	utils.TimeStats("Round", times["Round"])
	utils.TimeStats("LastRound", times["LastRound"])

	if err := blockpack.Check(ciP, ctx.EcdCI, ctx.DecCI, st, states); err != nil {
		t.Errorf("AES-128: %v", err)
	} else {
		fmt.Printf("\n=== OK: full AES-128, all %d blocks conform to the row-major AES oracle ===\n", len(states))
	}
}

// encRK packs a round key at the given level: the key stream is shared, so the 16 bytes are
// replicated over the whole batch and packed exactly like a state.
func encRK(t *testing.T, ctx *Context, rk []byte, blocks, level int) blockpack.Packed {
	t.Helper()
	repl := make([][16]byte, blocks)
	for bi := range repl {
		copy(repl[bi][:], rk)
	}
	p, err := blockpack.Encrypt(ctx.Sw.CiP, ctx.EcdCI, ctx.EncCI, repl, level)
	if err != nil {
		t.Fatalf("encrypt round key at level %d: %v", level, err)
	}
	return p
}

// cmp16 counts how many bytes differ between got and want and returns the first differing index.
func cmp16(got, want [16]byte) (wrong, first int) {
	first = -1
	for b := 0; b < 16; b++ {
		if got[b] != want[b] {
			if first < 0 {
				first = b
			}
			wrong++
		}
	}
	return
}

// report prints the slot + chain of st[0][0], the precision pooled over the 64 ciphertexts, and the
// whole batch against the oracle.
func report(ctx *Context, st blockpack.Packed, states [][16]byte, round int, step string) {
	ciP := ctx.Sw.CiP
	debug.DbgSlotCI(fmt.Sprintf("T%d %-11s st[0][0] =", round, step), st[0][0])
	debug.DbgChain(fmt.Sprintf("T%d %-11s chain    :", round, step), ctx.Sw.EvalCI, st[0][0])

	entries := make([]debug.PrecCt, 0, 64)
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			entries = append(entries, debug.PrecCt{
				Name: fmt.Sprintf("st[%2d][%d]", g, b),
				Want: blockpack.SlotVec(ciP, states, g, b),
				Ct:   st[g][b],
			})
		}
	}
	debug.PrecPoolCI(fmt.Sprintf("T%d %-11s prec (128 bits) :", round, step), entries...)

	got, bitErr, err := blockpack.Decrypt(ciP, ctx.EcdCI, ctx.DecCI, st)
	if err != nil {
		fmt.Printf("  [T%d] %-11s ORACLE: decrypt error: %v\n", round, step, err)
		return
	}
	wrong, firstBad := 0, -1
	for bi := range states {
		if got[bi] != states[bi] {
			wrong++
			if firstBad < 0 {
				firstBad = bi
			}
		}
	}
	if wrong == 0 {
		fmt.Printf("  [T%d] %-11s ORACLE: TRUE all %d blocks (worst bit err %.4f)\n", round, step, len(states), bitErr)
		return
	}
	w, f := cmp16(got[firstBad], states[firstBad])
	fmt.Printf("  [T%d] %-11s ORACLE: FALSE %d/%d blocks wrong (worst bit err %.4f; first bad block %d: %d/16 bytes, first byte %d)\n            got =%x\n            want=%x\n",
		round, step, wrong, len(states), bitErr, firstBad, w, f, got[firstBad], states[firstBad])
}
