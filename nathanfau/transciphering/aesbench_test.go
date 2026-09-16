package transciphering

//	go test ./nathanfau/transciphering/ -run '^TestAESBench$' -v -subbytes 2 -timeout 0
//	go test ./nathanfau/transciphering/ -run '^TestAESBench$' -v -timeout 0 -seed 42 -csv runs/aesbench.csv
//
// The -csv file gets ONE row per run, appended, so a file accumulates runs. Its columns are not
// those TestAES writes: give the two tests different files.

import (
	"fmt"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockpack"
	"github.com/tuneinsight/lattigo/v6/nathanfau/params2"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
)

// csvBenchRowHeader is what a bench run adds to the parameter columns: the two chronos, the memory,
// and the state as the oracle found it, once, at the end. kg_ms is the whole key generation, i.e.
// the keygen_ms and dft_ms columns plus the encryption of the key schedule. live_gb is the live heap
// once the keys are ready, peak_rss_gb the most the process ever held over the run. Sizes are in Go,
// 1 Go = 2^30 octets, like key_gb.
var csvBenchRowHeader = []string{
	"kg_ms", "transcipher_ms", "ms_per_block",
	"live_gb", "peak_rss_gb",
	"level", "prec_avg", "prec_min", "worst_err", "slots_pooled", "bit_err", "blocks_wrong",
}

var csvBenchHeader = append(append([]string{}, csvRunHeader...), csvBenchRowHeader...)

// TestAESBench is TestAES stripped down to what a timing should see: no trace, no per-operation
// oracle, the pipeline silenced, one line of progress per round. Two chronos, the key generation
// (FHE keys and the encrypted key schedule) and the transciphering itself (FirstRound -> 9 * Round
// -> LastRoundV1); the oracle runs once, after both have stopped.
func TestAESBench(t *testing.T) {
	const logN, k = 11, 4

	if testing.Short() {
		t.Skip("full AES: 10 rounds, several minutes")
	}

	cfg, sh := config(t), shape(t)
	seed := blockSeed()
	// The csv path is checked BEFORE the keys, so a bad path or a file with older columns fails in
	// a second instead of at the end of the run.
	if err := csvPrepare(*csvFlag, csvBenchHeader); err != nil {
		t.Fatalf("%v", err)
	}

	key := [16]byte{0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c}
	rk := aes.KeyExpansion(key[:])
	last := len(rk) - 1 // round 0 is the initial ARK, the last one has no MixColumns

	fmt.Printf(" AES-128 bench: SubBytes=V%d, XOR=%s, clean=%s, place=%s, extract=%s, STC=%s, logN=%d, random seed = %d \n",
		*sbVersion, cfg.Xor, cfg.Clean, cfg.Place, extractName(cfg), params2.FormatLogSTC(sh.LogSTC), logN, seed)

	// KeyGen: everything done once per key, i.e. the FHE keys and the key schedule encrypted at the
	// three levels the pipeline takes it at.
	tKG := time.Now()
	ctx, err := NewContextWith2(logN, k, cfg, sh)
	if err != nil {
		t.Fatalf("NewContext: %v", err)
	}
	ctx.Quiet = true
	ciP := ctx.Sw.CiP
	nBlocks := blockpack.Capacity(ciP)
	rkHE := make([]blockpack.Packed, len(rk))
	rkHE[0] = encRK(t, ctx, rk[0], nBlocks, ctx.InitLv)
	for r := 1; r < last; r++ {
		rkHE[r] = encRK(t, ctx, rk[r], nBlocks, ctx.ARKKeyLv)
	}
	rkHE[last] = encRK(t, ctx, rk[last], nBlocks, ctx.LastKeyLv)
	dKG := time.Since(tKG)
	fmt.Printf("KeyGen done in %s (bootstrapping %s, switcher %s, DFT matrices %s, key material %s)\n",
		dKG.Round(time.Millisecond), ctx.BtpKeyGen.Round(time.Millisecond), ctx.SwKeyGen.Round(time.Millisecond),
		ctx.EvalSetup.Round(time.Millisecond), utils.Bytes(uint64(ctx.BtpKeyBytes+ctx.SwKeyBytes)))
	// Between the chronos: the table collects.
	live := memTable("after KeyGen", ctx, packedItem("round keys", rkHE...))

	blocks := blockpack.RandomBlocks(ciP, rand.New(rand.NewSource(seed)))

	// Transciphering: the blocks arrive in the clear, FirstRound encodes them itself. One line per
	// round, and nothing else: a round costs minutes, a print costs microseconds.
	tTC := time.Now()
	st, err := ctx.FirstRound(blocks, rkHE[0])
	if err != nil {
		t.Fatalf("FirstRound: %v", err)
	}
	for r := 1; r <= last; r++ {
		tRound := time.Now()
		if r < last {
			st, err = ctx.Round(st, rkHE[r], *sbVersion, nil)
		} else {
			st, err = ctx.LastRoundV1(st, rkHE[r], *sbVersion, nil)
		}
		if err != nil {
			t.Fatalf("Round T%d: %v", r, err)
		}
		fmt.Printf("  round %2d/%d done in %s (elapsed %s, level %d) | %s\n",
			r, last, time.Since(tRound).Round(time.Millisecond), time.Since(tTC).Round(time.Millisecond), st[0][0].Level(), utils.Mem())
	}
	dTC := time.Since(tTC)

	memTable("after the rounds", ctx, packedItem("round keys", rkHE...), packedItem("state", st))
	peak := utils.Mem().PeakRSS
	fmt.Printf("KeyGen         : %s\n", dKG.Round(time.Millisecond))
	fmt.Printf("Transciphering : %s  (%d blocks, %s per block)\n",
		dTC.Round(time.Millisecond), nBlocks, (dTC / time.Duration(nBlocks)).Round(time.Microsecond))
	fmt.Printf("Memory         : live %s after KeyGen, peak RSS %s over the run\n", utils.Bytes(live), utils.Bytes(peak))

	// The one oracle call: the whole cipher in the clear, then the state decrypted once.
	want := aesRM(blocks, rk)
	prec, got, err := measure(ctx, st, want)
	if err != nil {
		t.Fatalf("final measure: %v", err)
	}
	fmt.Printf("Precision      : avg %.2f bits, worst %.2f bits over %d slots, worst bit err %.4f\n",
		prec.AvgPrec, negLog2(prec.WorstBit), prec.Slots, prec.BitErr)

	rec, err := newAESBenchCSV(*csvFlag, ctx, cfg, logN, k, nBlocks, last-1, seed)
	if err != nil {
		t.Errorf("%v", err)
	}
	rec.addBench(dKG, dTC, nBlocks, live, peak, prec)
	if err := rec.write(); err != nil {
		t.Errorf("%v", err)
	}

	if prec.Wrong != 0 {
		bad := firstDiff(got, want)
		w, f := cmp16(got[bad], want[bad])
		t.Errorf("AES-128: %d/%d blocks wrong (first bad block %d: %d/16 bytes, first byte %d)\n  got =%x\n  want=%x",
			prec.Wrong, nBlocks, bad, w, f, got[bad], want[bad])
		return
	}
	fmt.Printf("=== OK: full AES-128, all %d blocks conform to the row-major AES oracle ===\n", nBlocks)
}

// aesRM is the cleartext cipher the pipeline computes, on blockpack's row-major layout: the
// sequence TestAES advances step by step, run in one go. The round keys are XORed as they come out
// of KeyExpansion, exactly as encRK packs them.
func aesRM(blocks [][16]byte, rk [][]byte) [][16]byte {
	last := len(rk) - 1
	out := make([][16]byte, len(blocks))
	for bi, s := range blocks {
		aes.AddRoundKey(s[:], rk[0])
		for r := 1; r <= last; r++ {
			aes.SubBytes(s[:])
			s = aes.ShiftRowsRM(s)
			if r < last {
				s = aes.MixColumnsRM(s)
			}
			aes.AddRoundKey(s[:], rk[r])
		}
		out[bi] = s
	}
	return out
}

// newAESBenchCSV is newAESCSV on the bench columns: the same parameter columns, one row per run
// instead of one per operation.
func newAESBenchCSV(path string, ctx *Context, cfg Config, logN, k, blocks, rounds int, seed int64) (*aesCSV, error) {
	if path == "" {
		return nil, nil
	}
	run, err := csvRunValues(ctx, cfg, logN, k, blocks, rounds, seed)
	if err != nil {
		return nil, err
	}
	return &aesCSV{path: path, header: csvBenchHeader, run: run}, nil
}

// addBench records the single row of a bench run. prec.WorstBit is a distance, so prec_min is its
// -log2: the MINIMUM precision over the pool, the number that decides whether a bit is about to
// flip.
func (c *aesCSV) addBench(kg, tc time.Duration, blocks int, live, peak uint64, prec precSummary) {
	if c == nil {
		return
	}
	gb := func(b uint64) string { return ftoa(float64(b)/(1<<30), 3) }
	row := append([]string{}, c.run...)
	row = append(row,
		ftoa(millis(kg), 3), ftoa(millis(tc), 3), ftoa(millis(tc)/float64(blocks), 3),
		gb(live), gb(peak),
		itoa(prec.Level), ftoa(prec.AvgPrec, 3), ftoa(negLog2(prec.WorstBit), 3), ftoa(prec.WorstBit, 6),
		itoa(prec.Slots), ftoa(prec.BitErr, 6), itoa(prec.Wrong))
	c.rows = append(c.rows, row)
}
