package transciphering

//	go test ./nathanfau/transciphering/ -run '^TestTransciphering$' -v -subbytes 2 -rounds 1 -timeout 0
//	go test ./nathanfau/transciphering/ -run '^TestTransciphering$' -v -subbytes 2 -rounds 4 -timeout 0 -csv runs/aes.csv
//	go test ./nathanfau/transciphering/ -run '^TestAES$' -v -subbytes 2 -timeout 0
//	go test ./nathanfau/transciphering/ -run '^TestAES$' -v -subbytes 2 -timeout 0 -logn 12 -stc 60 -logqi 40

import (
	"flag"
	"fmt"
	"math/rand"
	"strconv"

	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockpack"
	"github.com/tuneinsight/lattigo/v6/nathanfau/cleaning"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/nathanfau/params2"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
)

var (
	nRoundsFlag  = flag.Int("rounds", 1, "number of AES middle rounds to run")
	sbVersion    = flag.Int("subbytes", 2, "SubBytes version (1..3), by decreasing cost: 247, 98, 69 relin per byte")
	xorFlag      = flag.String("xor", "nosq", `XOR circuit used throughout: "nosq" for x+y-2xy, "sq" for (x-y)^2`)
	cleanFlag    = flag.String("clean", "cleaning", `cleaning polynomial: "cleaning" (2 levels), "smoother" or "verysmoother" (3 levels, one prime more)`)
	placeFlag    = flag.String("place", "after", `where the round cleaning lands: "after" (AddRoundKey then Cleaning), "both" (clean state and key, then XOR) or "one" (clean the state only, then XOR)`)
	xtractFlag   = flag.Bool("cleanextract", false, "run the refresh on BitExtractClean (interpolation and cleaning fused, error quadratic in Algo1's output) instead of BitExtract (half spectrum, error linear); costs one prime more")
	seedFlag     = flag.Int64("seed", 0, "seed of the block draw; 0 draws one from the clock")
	stcFlag      = flag.String("stc", "2x30", `SlotsToCoeffs primes, one level each: "2x30" or "30,30"; "60" is the chain's historical single prime (one level, dense matrix, far too heavy at large logN)`)
	logNFlag     = flag.Int("logn", 11, "ring degree of the pipeline")
	logQiFlag    = flag.Int("logqi", params2.LogScale, "size in bits of the chain's primes, every Q prime but q0 and the SlotsToCoeffs ones, which is also the scale the pipeline runs at; q0 follows it, P and -stc do not")
	sbExactFlag  = flag.Bool("sbexact", false, "run SubBytes on the exactly aligned S-box (aes.SubByteExact); temporary. Implies -cleanfixed")
	refreshCanon = flag.Bool("refreshcanon", false, "land the refresh's bit extraction on the canonical scale instead of relabelling its output as canonical")
	zoneFlag     = flag.Int("zone", 0, "size of the primes of the AES circuit and the bit extraction, and of the scale there (params2.Shape.LogZone); 0 or -logqi = no zone, the pipeline as it always ran")
	cleanFixed   = flag.Bool("cleanfixed", false, "land the round cleaning on the default scale instead of its input's, so the 8 bits of a byte leave it on one scale; -sbexact implies it")
)

// sboxName is the SubBytes variant the flags select, as the headers and the csv write it.
func sboxName() string {
	if *sbExactFlag {
		return fmt.Sprintf("%dexact", *sbVersion)
	}
	return strconv.Itoa(*sbVersion)
}

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
	return Config{Xor: xk, Clean: ck, Place: pk, CleanExtract: *xtractFlag}
}

// shape parses -stc, -logqi and -zone, also before any keygen. Every level the tests use is read off the
// context, which follows the extra levels a longer SlotsToCoeffs adds.
func shape(t *testing.T) params2.Shape {
	t.Helper()
	logSTC, err := params2.ParseLogSTC(*stcFlag)
	if err != nil {
		t.Fatalf("-stc: %v", err)
	}
	if *logQiFlag < 1 {
		t.Fatalf("-logqi %d: want a prime size in bits", *logQiFlag)
	}
	return params2.Shape{LogSTC: logSTC, LogQi: *logQiFlag, LogZone: *zoneFlag}
}

// chainName is what the headers print of the chain the flags selected.
func chainName(sh params2.Shape) string {
	zone := "none"
	if sh.HasZone() {
		zone = fmt.Sprintf("%d bits", sh.LogZone)
	}
	return fmt.Sprintf("logN=%d, q_i=%d, STC=%s, zone=%s", *logNFlag, sh.LogQi, params2.FormatLogSTC(sh.LogSTC), zone)
}

// extractName is what the trace calls the extraction the config selected.
func extractName(cfg Config) string {
	if cfg.CleanExtract {
		return "clean (interp+cleaning, k+1 lv)"
	}
	return "half (demi-spectre, k lv)"
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
	const k = 4
	logN := *logNFlag

	n := *nRoundsFlag
	if n < 1 {
		n = 1
	}
	cfg, sh := config(t), shape(t)
	fmt.Printf(" Transciphering: rounds=%d, SubBytes=V%s, XOR=%s, clean=%s, place=%s, extract=%s, %s \n",
		n, sboxName(), cfg.Xor, cfg.Clean, cfg.Place, extractName(cfg), chainName(sh))

	ctx, err := NewContextWith2(logN, k, cfg, sh)
	if err != nil {
		t.Fatalf("NewContext: %v", err)
	}
	ctx.SBoxExact = *sbExactFlag
	ctx.CleanFixedScale = *cleanFixed
	ctx.RefreshCanon = *refreshCanon
	debug.DbgParams("TranscipheringParams", ctx.Params)
	ciP := ctx.Sw.CiP

	// FIPS-197 key (shared key stream). Each slot carries a DISTINCT random block (real batching);
	// states[bi] tracks block bi's cleartext state, and entry to round 1 is plaintext XOR rk0.
	key := [16]byte{0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c}
	rk := aes.KeyExpansion(key[:])

	seed := blockSeed()
	states := blockpack.RandomBlocks(ciP, rand.New(rand.NewSource(seed)))
	fmt.Printf("random seed = %d  (%d blocks)\n", seed, len(states))

	rec, err := newAESCSV(*csvFlag, ctx, cfg, logN, k, len(states), n, seed)
	if err != nil {
		t.Fatalf("%v", err)
	}
	for bi := range states {
		aes.AddRoundKey(states[bi][:], rk[0])
	}

	st, err := blockpack.EncryptAt(ciP, ctx.EcdCI, ctx.EncCI, states, ctx.SubBytesLv, ctx.Canon)
	if err != nil {
		t.Fatalf("blockpack.Encrypt state: %v", err)
	}

	rkHE := make([]blockpack.Packed, n+1)
	for r := 1; r <= n; r++ {
		rkHE[r] = encRK(t, ctx, rk[r], len(states), ctx.ARKKeyLv)
	}
	memTable("before the rounds", ctx, packedItem("round keys", rkHE...), packedItem("state", st))

	fmt.Println("================ input (entry to round 1) ================")
	rec.add(0, "input", 0, report(ctx, st, states, 0, "input"))

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
			rec.add(r, name, dur, report(ctx, st, states, r, name))
		}

		if st, err = ctx.Round(st, rkHE[r], *sbVersion, after); err != nil {
			t.Fatalf("Round T%d: %v", r, err)
		}

		dRound := time.Since(tRound)
		times["Round"] = append(times["Round"], dRound)
		fmt.Printf("  [T%d] round time (incl. oracle traces): %s\n", r, dRound.Round(time.Millisecond))
		memLine(fmt.Sprintf("after T%d", r))
	}
	memTable("after the rounds", ctx, packedItem("round keys", rkHE...), packedItem("state", st))

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

	if err := rec.write(); err != nil {
		t.Errorf("%v", err)
	}
}

// TestAES runs the whole cipher on a full batch of DISTINCT random blocks:
// FirstRound -> 9 * Round -> LastRoundV1, with the oracle checked after EVERY operation.
func TestAES(t *testing.T) {
	const k = 4
	logN := *logNFlag

	if testing.Short() {
		t.Skip("full AES: 10 rounds, several minutes")
	}

	cfg, sh := config(t), shape(t)

	ctx, err := NewContextWith2(logN, k, cfg, sh)
	if err != nil {
		t.Fatalf("NewContext: %v", err)
	}
	ctx.SBoxExact = *sbExactFlag
	ctx.CleanFixedScale = *cleanFixed
	ctx.RefreshCanon = *refreshCanon
	debug.DbgParams("TranscipheringParams", ctx.Params)
	ciP := ctx.Sw.CiP

	key := [16]byte{0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c}
	rk := aes.KeyExpansion(key[:])
	nMiddle := len(rk) - 2 // round 0 is the initial ARK, the last one has no MixColumns

	seed := blockSeed()
	blocks := blockpack.RandomBlocks(ciP, rand.New(rand.NewSource(seed)))
	fmt.Printf(" AES-128: %d middle rounds, SubBytes=V%s, XOR=%s, clean=%s, place=%s, extract=%s, %s, %d blocks, random seed = %d \n",
		nMiddle, sboxName(), cfg.Xor, cfg.Clean, cfg.Place, extractName(cfg), chainName(sh), len(blocks), seed)

	rec, err := newAESCSV(*csvFlag, ctx, cfg, logN, k, len(blocks), nMiddle, seed)
	if err != nil {
		t.Fatalf("%v", err)
	}

	// The blocks stay in the CLEAR: that is what the client sends and what FirstRound takes.
	states := make([][16]byte, len(blocks))
	copy(states, blocks)

	// Three key levels: rk0 lands before SubBytes, the middle keys after MixColumns, the last one
	// straight out of the refresh since no MixColumns eats into it.
	rk0 := encRK(t, ctx, rk[0], len(blocks), ctx.InitLv)
	rkHE := make([]blockpack.Packed, len(rk))
	for r := 1; r <= nMiddle; r++ {
		rkHE[r] = encRK(t, ctx, rk[r], len(blocks), ctx.ARKKeyLv)
	}
	rkHE[len(rk)-1] = encRK(t, ctx, rk[len(rk)-1], len(blocks), ctx.LastKeyLv)
	memTable("before the rounds", ctx, packedItem("round keys", append(rkHE, rk0)...))

	times := map[string][]time.Duration{}
	tGlobal := time.Now()

	fmt.Println("\n================ Round 0 (initial ARK) ================")
	t0 := time.Now()
	st, err := ctx.FirstRound(blocks, rk0)
	if err != nil {
		t.Fatalf("FirstRound: %v", err)
	}
	dFirst := time.Since(t0)
	times["FirstRound"] = append(times["FirstRound"], dFirst)
	for bi := range states {
		aes.AddRoundKey(states[bi][:], rk[0])
	}
	if l := st[0][0].Level(); l != ctx.SubBytesLv {
		t.Errorf("FirstRound left the state at level %d, want the SubBytes level %d", l, ctx.SubBytesLv)
	}
	rec.add(0, "FirstRound", dFirst, report(ctx, st, states, 0, "FirstRound"))
	memLine("after T0")

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
		rec.add(round, name, dur, report(ctx, st, states, round, name))
	}

	for r := 1; r <= nMiddle; r++ {
		fmt.Printf("\n================ Round %d/%d ================\n", r, len(rk)-1)
		round, rkNow = r, rk[r]
		tRound := time.Now()
		if st, err = ctx.Round(st, rkHE[r], *sbVersion, after); err != nil {
			t.Fatalf("Round T%d: %v", r, err)
		}
		times["Round"] = append(times["Round"], time.Since(tRound))
		memLine(fmt.Sprintf("after T%d", r))
	}

	last := len(rk) - 1
	fmt.Printf("\n================ Round %d/%d (last, no MixColumns) ================\n", last, last)
	round, rkNow = last, rk[last]
	tRound := time.Now()
	if st, err = ctx.LastRoundV1(st, rkHE[last], *sbVersion, after); err != nil {
		t.Fatalf("LastRoundV1: %v", err)
	}
	times["LastRound"] = append(times["LastRound"], time.Since(tRound))
	memLine(fmt.Sprintf("after T%d", last))
	memTable("after the rounds", ctx, packedItem("round keys", append(rkHE, rk0)...), packedItem("state", st))

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

	if err := rec.write(); err != nil {
		t.Errorf("%v", err)
	}
}

// encRK packs a round key at the given level: the key stream is shared, so the 16 bytes are
// replicated over the whole batch and packed exactly like a state, at the state's scale.
func encRK(t *testing.T, ctx *Context, rk []byte, blocks, level int) blockpack.Packed {
	t.Helper()
	repl := make([][16]byte, blocks)
	for bi := range repl {
		copy(repl[bi][:], rk)
	}
	p, err := blockpack.EncryptAt(ctx.Sw.CiP, ctx.EcdCI, ctx.EncCI, repl, level, ctx.Canon)
	if err != nil {
		t.Fatalf("encrypt round key at level %d: %v", level, err)
	}
	return p
}

// packedItem weighs packed states or round keys as one line of memTable; unset ones are skipped.
func packedItem(name string, ps ...blockpack.Packed) MemItem {
	it := MemItem{Name: name}
	for _, p := range ps {
		for _, ct := range p.Cts() {
			if ct != nil {
				it.Count++
			}
		}
		it.Bytes += p.BinarySize()
	}
	return it
}

// memTable prints what the run keeps alive, object by object, against the live heap, and returns
// that heap. The gap is what no object reports: evaluator buffers, encoders, ring tables, Go
// overhead. Two tables of one run that differ in live heap beyond the listed objects point at
// something the rounds retain. It collects, so keep it out of timed sections.
func memTable(title string, ctx *Context, extra ...MemItem) uint64 {
	items := append(ctx.MemBreakdown(), extra...)
	live := utils.LiveHeap()
	share := func(b uint64) string { return fmt.Sprintf("%5.1f%%", 100*float64(b)/float64(live)) }

	fmt.Printf("---------------- memory %s ----------------\n", title)
	var listed uint64
	for _, it := range items {
		listed += uint64(it.Bytes)
		fmt.Printf("  %-26s %5d  %12s  %s\n", it.Name, it.Count, utils.Bytes(uint64(it.Bytes)), share(uint64(it.Bytes)))
	}
	gap := "  (listed > live)"
	if live >= listed {
		gap = share(live - listed)
	}
	fmt.Printf("  %-26s %5s  %12s  %s\n", "not listed (buffers...)", "", utils.BytesDelta(live, listed), gap)
	fmt.Printf("  %-26s %5s  %12s\n", "live heap", "", utils.Bytes(live))
	fmt.Printf("  %s | %s\n", utils.GCSettings(), utils.Mem())
	return live
}

// memLine is one line of memory state, cheap enough for a timed section: it does not collect.
func memLine(label string) {
	fmt.Printf("  [mem] %-12s %s\n", label, utils.Mem())
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
// precSummary is what report measured, so a caller can write it out instead of only reading it.
type precSummary struct {
	AvgPrec  float64 // mean precision in bits, pooled over the 64 ciphertexts and all their slots
	WorstBit float64 // largest |value - reference| over the same pool, i.e. the MINIMUM precision
	Slots    int     // how many values that pool holds
	Level    int
	Wrong    int     // blocks the oracle disagrees with
	BitErr   float64 // worst distance to a clean bit, from blockpack.Decrypt
}

// measure is what a row records, and nothing else: the precision pooled over the 64 ciphertexts and
// the oracle verdict, with no printing, so a timed run can call it once the chrono has stopped. It
// also hands back the decrypted blocks, so a caller can name the first one that went wrong without
// decrypting a second time.
func measure(ctx *Context, st blockpack.Packed, want [][16]byte) (precSummary, [][16]byte, error) {
	ciP := ctx.Sw.CiP
	sum := precSummary{Level: st[0][0].Level()}

	stats := make([]utils.BitStats, 0, 64)
	for g := 0; g < 8; g++ {
		for b := 0; b < 8; b++ {
			s, err := utils.BitDistanceCt(ctx.EcdCI, ctx.DecCI, st[g][b], blockpack.SlotVec(ciP, want, g, b), 0)
			if err != nil {
				return sum, nil, fmt.Errorf("prec st[%d][%d]: %w", g, b, err)
			}
			stats = append(stats, s)
		}
	}
	agg, _ := utils.WorstOf(stats)
	sum.AvgPrec, sum.WorstBit, sum.Slots = agg.AvgPrec, agg.Worst, agg.Slots

	got, bitErr, err := blockpack.Decrypt(ciP, ctx.EcdCI, ctx.DecCI, st)
	if err != nil {
		return sum, nil, fmt.Errorf("decrypt: %w", err)
	}
	sum.BitErr = bitErr
	for bi := range want {
		if got[bi] != want[bi] {
			sum.Wrong++
		}
	}
	return sum, got, nil
}

// firstDiff is the index of the first block the oracle disagrees with, -1 when there is none.
func firstDiff(got, want [][16]byte) int {
	for bi := range want {
		if got[bi] != want[bi] {
			return bi
		}
	}
	return -1
}

func report(ctx *Context, st blockpack.Packed, states [][16]byte, round int, step string) precSummary {
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

	// The same pool again, kept this time: one avg and one worst over ALL 64 ciphertexts and all
	// their slots, not per byte, plus the batch against the oracle.
	sum, got, err := measure(ctx, st, states)
	if err != nil {
		fmt.Printf("  [T%d] %-11s %v\n", round, step, err)
		return sum
	}
	if sum.Wrong == 0 {
		fmt.Printf("  [T%d] %-11s ORACLE: TRUE all %d blocks (worst bit err %.4f)\n", round, step, len(states), sum.BitErr)
		return sum
	}
	firstBad := firstDiff(got, states)
	w, f := cmp16(got[firstBad], states[firstBad])
	fmt.Printf("  [T%d] %-11s ORACLE: FALSE %d/%d blocks wrong (worst bit err %.4f; first bad block %d: %d/16 bytes, first byte %d)\n            got =%x\n            want=%x\n",
		round, step, sum.Wrong, len(states), sum.BitErr, firstBad, w, f, got[firstBad], states[firstBad])
	return sum
}
