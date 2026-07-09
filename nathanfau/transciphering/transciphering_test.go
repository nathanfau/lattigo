package transciphering

// Run examples:
//   go test ./nathanfau/transciphering/ -run TestTransciphering -v
//   go test ./nathanfau/transciphering/ -run TestTransciphering -v -rounds 2
//   go test ./nathanfau/transciphering/ -run TestTransciphering -v -subbytes 3 -mixcolumns 1

import (
	"flag"
	"fmt"
	"sort"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
)

var (
	nRoundsFlag = flag.Int("rounds", 1, "number of AES middle rounds to run")
	sbVersion   = flag.Int("subbytes", 2, "SubBytes version (1..4)")
	mcVersion   = flag.Int("mixcolumns", 2, "MixColumns version (1..2)")
	noTrickFlag = flag.Bool("notrick", false, "use the paper-faithful AES evaluator (no squaring/scale-chain trick)")
)

func to16(b []byte) [16]byte {
	var out [16]byte
	copy(out[:], b)
	return out
}

// timestats logs min/max/mean/median over a list of per-round durations for one step.
func timestats(t *testing.T, name string, ds []time.Duration) {
	t.Helper()
	if len(ds) == 0 {
		t.Logf("  %-13s (no measure)", name)
		return
	}
	cp := append([]time.Duration(nil), ds...)
	sort.Slice(cp, func(i, j int) bool { return cp[i] < cp[j] })
	var sum time.Duration
	for _, d := range cp {
		sum += d
	}
	mean := sum / time.Duration(len(cp))
	med := cp[len(cp)/2]
	t.Logf("  %-13s min=%s max=%s mean=%s med=%s (n=%d)", name,
		cp[0].Round(time.Millisecond), cp[len(cp)-1].Round(time.Millisecond),
		mean.Round(time.Millisecond), med.Round(time.Millisecond), len(cp))
}

// cmp16 counts wrong bytes and returns the first differing index (-1 if equal).
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

// TestTransciphering runs n AES middle rounds (SubBytes -> refresh(Algo1) -> ShiftRows ->
// MixColumns -> AddRoundKey -> Cleaning) on an encrypted state, with slot/chain/precision
// traces and an AES-oracle comparison after every operation (intra-round and between rounds).
//
// The number of rounds and the SubBytes / MixColumns versions are set by -rounds, -subbytes
// and -mixcolumns.
func TestTransciphering(t *testing.T) {
	const logN, k = 11, 4

	n := *nRoundsFlag
	if n < 1 {
		n = 1
	}
	fmt.Printf("=== Transciphering: rounds=%d, SubBytes=V%d, MixColumns=V%d, noTrick=%v ===\n", n, *sbVersion, *mcVersion, *noTrickFlag)

	ctx, err := NewContext(logN, k, *noTrickFlag, *mcVersion)
	if err != nil {
		t.Fatalf("NewContext: %v", err)
	}

	// FIPS-197 key / plaintext, used to drive both the HE circuit and the cleartext oracle.
	key := [16]byte{0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c}
	pt := [16]byte{0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34}
	rk := aes.KeyExpansion(key[:])

	// Oracle state = entry to round 1 (pt XOR rk0). Rounds 1..n are full middle rounds.
	want := pt
	aes.AddRoundKey(want[:], rk[0])

	// Input state at the loop level [42, 60, 38, 38, 38, 38] = level 5: 3 levels for SubBytes,
	// 1 for the CIToStandard conversion, 1 (the 60) for SlotsToCoeffs, before the refresh recharges.
	const inputLevel = 5
	keyLevel := ctx.Sw.CiP.MaxLevel()
	st, err := aes.EncStateRepl(ctx.Sw.CiP, ctx.EcdCI, ctx.EncCI, want, inputLevel)
	if err != nil {
		t.Fatalf("EncStateRepl state: %v", err)
	}
	rkHE := make([]aes.StateHE, n+1)
	for r := 1; r <= n; r++ {
		if rkHE[r], err = aes.EncStateRepl(ctx.Sw.CiP, ctx.EcdCI, ctx.EncCI, to16(rk[r]), keyLevel); err != nil {
			t.Fatalf("EncStateRepl rk%d: %v", r, err)
		}
	}

	// dec0 decrypts the state and returns slot 0 (representative of the replicated batch).
	dec0 := func(s aes.StateHE) [16]byte {
		outs, derr := aes.DecStateBatch(ctx.Sw.CiP, ctx.EcdCI, ctx.DecCI, s)
		if derr != nil {
			t.Fatalf("DecStateBatch: %v", derr)
		}
		return outs[0]
	}

	// wantBitSlice returns the (replicated) expected value of one state bit for the precision trace.
	slots := ctx.Sw.CiP.MaxSlots()
	wantBitSlice := func(byteIdx, bit int) []float64 {
		v := float64((want[byteIdx] >> uint(bit)) & 1)
		s := make([]float64, slots)
		for i := range s {
			s[i] = v
		}
		return s
	}

	// report dumps slot + chain + precision on two representative bits and compares the whole
	// decrypted state against the current oracle state.
	report := func(round int, step string) {
		debug.DbgSlotCI(fmt.Sprintf("T%d %-11s st[0][0] =", round, step), st[0][0])
		debug.DbgChain(fmt.Sprintf("T%d %-11s chain    :", round, step), ctx.Sw.EvalCI, st[0][0])
		debug.PrintPrecCI(fmt.Sprintf("T%d %-11s prec[0][0] :", round, step), wantBitSlice(0, 0), st[0][0])
		debug.DbgSlotCI(fmt.Sprintf("T%d %-11s st[7][4] =", round, step), st[7][4])
		debug.PrintPrecCI(fmt.Sprintf("T%d %-11s prec[7][4] :", round, step), wantBitSlice(7, 4), st[7][4])

		got := dec0(st)
		if w, f := cmp16(got, want); w > 0 {
			fmt.Printf("  [T%d] %-11s ORACLE: FALSE %2d/16 bytes wrong (first = byte %d)\n            got =%x\n            want=%x\n", round, step, w, f, got, want)
		} else {
			fmt.Printf("  [T%d] %-11s ORACLE: TRUE 16/16\n", round, step)
		}
	}

	fmt.Println("================ input (entry to round 1) ================")
	report(0, "input")

	// Timers: per-step durations accumulated across rounds (+ per-round total), plus a global one.
	times := map[string][]time.Duration{}
	tGlobal := time.Now()

	for r := 1; r <= n; r++ {
		fmt.Printf("\n================ Round %d/%d ================\n", r, n)
		tRound := time.Now()

		// SubBytes
		fmt.Printf("---- T%d SubBytes (V%d) ----\n", r, *sbVersion)
		t0 := time.Now()
		if st, err = ctx.AE.SubBytes(st, *sbVersion); err != nil {
			t.Fatalf("SubBytes T%d: %v", r, err)
		}
		times["SubBytes"] = append(times["SubBytes"], time.Since(t0))
		aes.SubBytes(want[:])
		report(r, "SubBytes")

		// Refresh (Algo1). The bootstrap preserves the bit values, so the oracle is unchanged.
		fmt.Printf("---- T%d Refresh (Algo1) ----\n", r)
		t0 = time.Now()
		if st, err = ctx.RefreshState(st); err != nil {
			t.Fatalf("RefreshState T%d: %v", r, err)
		}
		times["Refresh"] = append(times["Refresh"], time.Since(t0))
		report(r, "Refresh")

		// ShiftRows
		fmt.Printf("---- T%d ShiftRows ----\n", r)
		t0 = time.Now()
		st = aes.ShiftRowsHE(st)
		times["ShiftRows"] = append(times["ShiftRows"], time.Since(t0))
		aes.ShiftRows(want[:])
		report(r, "ShiftRows")

		// MixColumns
		fmt.Printf("---- T%d MixColumns (V%d) ----\n", r, *mcVersion)
		t0 = time.Now()
		if st, err = ctx.AE.MixColumns(st, *mcVersion); err != nil {
			t.Fatalf("MixColumns T%d: %v", r, err)
		}
		times["MixColumns"] = append(times["MixColumns"], time.Since(t0))
		aes.MixColumns(want[:])
		report(r, "MixColumns")

		// AddRoundKey. The fresh round key is at MaxLevel while the state is deep, so it must be
		// canonicalized (DropLevel to the state level) first: otherwise trick-mode xor aligns by
		// squaring the key ~(MaxLevel-stateLevel) times and blows up its noise.
		fmt.Printf("---- T%d AddRoundKey ----\n", r)
		t0 = time.Now()
		if st, err = ctx.AddRoundKeyCanon(st, rkHE[r]); err != nil {
			t.Fatalf("AddRoundKey T%d: %v", r, err)
		}
		times["AddRoundKey"] = append(times["AddRoundKey"], time.Since(t0))
		aes.AddRoundKey(want[:], rk[r])
		report(r, "AddRoundKey")

		// Cleaning. Snaps bits to 0/1, oracle unchanged.
		fmt.Printf("---- T%d Cleaning ----\n", r)
		t0 = time.Now()
		if st, err = ctx.CleanState(st); err != nil {
			t.Fatalf("Cleaning T%d: %v", r, err)
		}
		times["Cleaning"] = append(times["Cleaning"], time.Since(t0))
		report(r, "Cleaning")

		dRound := time.Since(tRound)
		times["Round"] = append(times["Round"], dRound)
		fmt.Printf("  [T%d] round time (excl. oracle traces): %s\n", r, dRound.Round(time.Millisecond))
	}

	fmt.Printf("\nTOTAL (%d rounds, incl. oracle traces): %s\n", n, time.Since(tGlobal).Round(time.Millisecond))
	fmt.Println("================ timing stats (per round) ================")
	for _, name := range []string{"SubBytes", "Refresh", "ShiftRows", "MixColumns", "AddRoundKey", "Cleaning", "Round"} {
		timestats(t, name, times[name])
	}

	got := dec0(st)
	if got != want {
		t.Errorf("AES HE (%d rounds): got=%x want=%x", n, got, want)
	} else {
		fmt.Printf("\n=== OK: %d middle rounds conform to the AES oracle (16/16) ===\n", n)
	}
}
