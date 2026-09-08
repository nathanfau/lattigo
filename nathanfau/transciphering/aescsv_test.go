package transciphering

import (
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"slices"
	"strconv"
	"strings"
	"time"
)

// One row per operation of a TestAES run, APPENDED to -csv.
//
// Every row repeats the run's parameters, and carries the run's start time. That is deliberate:
// successive runs can share one file, a comparison is a group-by on run_ts or on whichever
// parameter changed, and no row is ever ambiguous about which run produced it.
//
//	go test ./nathanfau/transciphering/ -run '^TestAES$' -v -timeout 0 -seed 42 -csv runs/aes.csv
//	go test ./nathanfau/transciphering/ -run '^TestAES$' -v -timeout 0 -seed 42 -cleanextract -csv runs/aes.csv

var csvFlag = flag.String("csv", "", "append one row per operation of TestAES to this file: parameters, timing, and the precision pooled over all 64 ciphertexts")

var csvHeader = []string{
	// which run
	"run_ts",
	// what was run
	"logn", "k", "slots", "blocks", "seed", "rounds",
	"subbytes", "xor", "clean", "cleandepth", "place", "extract", "extractlv",
	"primes", "logscale", "refresh_lv", "ark_lv",
	// where in the run
	"seq", "round", "step",
	// what it cost
	"ms",
	// what came out
	"level", "prec_avg", "prec_min", "worst_err", "slots_pooled", "bit_err", "blocks_wrong",
}

// aesCSV accumulates the rows and writes them once, at the end: a run that dies mid-way still
// leaves nothing half-written to misread later.
type aesCSV struct {
	path string
	run  []string // the parameter columns, identical on every row
	rows [][]string
	seq  int
}

// newAESCSV also checks, before anything is computed, that the file can be appended to: finding out
// at the END of a ten-minute run that its header belongs to an older set of columns would be a poor
// way to learn it.
func newAESCSV(path string, ctx *Context, cfg Config, logN, k, blocks, rounds int, seed int64) (*aesCSV, error) {
	if path == "" {
		return nil, nil
	}
	if dir := filepath.Dir(path); dir != "." {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			return nil, fmt.Errorf("csv: %w", err)
		}
	}
	if _, err := csvIsFresh(path); err != nil {
		return nil, err
	}
	extractLv := cfg.ExtractLevels(k)
	return &aesCSV{
		path: path,
		run: []string{
			time.Now().Format(time.RFC3339),
			itoa(logN), itoa(k), itoa(ctx.Sw.CiP.MaxSlots()), itoa(blocks), strconv.FormatInt(seed, 10), itoa(rounds),
			itoa(*sbVersion), cfg.Xor.String(), cfg.Clean.String(), itoa(cfg.Clean.Depth()),
			cfg.Place.String(), extractName(cfg), itoa(extractLv),
			itoa(ctx.Params.MaxLevel() + 1), itoa(ctx.Params.LogDefaultScale()),
			itoa(ctx.RefreshLv), itoa(ctx.ARKLv),
		},
	}, nil
}

// add records one operation. prec.WorstBit is a distance, so prec_min is its -log2: the MINIMUM
// precision over the pool, the number that decides whether a bit is about to flip.
func (c *aesCSV) add(round int, step string, dur time.Duration, prec precSummary) {
	if c == nil {
		return
	}
	c.seq++
	row := append([]string{}, c.run...)
	row = append(row,
		itoa(c.seq), itoa(round), step,
		ftoa(float64(dur.Microseconds())/1000, 3),
		itoa(prec.Level), ftoa(prec.AvgPrec, 3), ftoa(negLog2(prec.WorstBit), 3), ftoa(prec.WorstBit, 6),
		itoa(prec.Slots), ftoa(prec.BitErr, 6), itoa(prec.Wrong))
	c.rows = append(c.rows, row)
}

// write APPENDS to the file, so successive runs accumulate: point every run at the same -csv and
// the comparison is a group-by on run_ts, or on whichever parameter column changed. The header goes
// in only when the file is new, and an existing file whose header is not this one is refused rather
// than corrupted -- that happens when the columns here change under an older file.
func (c *aesCSV) write() error {
	if c == nil {
		return nil
	}

	for i, r := range c.rows {
		if len(r) != len(csvHeader) {
			return fmt.Errorf("csv: row %d has %d fields for %d columns", i, len(r), len(csvHeader))
		}
	}

	fresh, err := csvIsFresh(c.path)
	if err != nil {
		return err
	}

	f, err := os.OpenFile(c.path, os.O_CREATE|os.O_WRONLY|os.O_APPEND, 0o644)
	if err != nil {
		return fmt.Errorf("csv: %w", err)
	}
	defer f.Close()

	w := csv.NewWriter(f)
	if fresh {
		if err := w.Write(csvHeader); err != nil {
			return fmt.Errorf("csv header: %w", err)
		}
	}
	if err := w.WriteAll(c.rows); err != nil {
		return fmt.Errorf("csv rows: %w", err)
	}
	w.Flush()
	if err := w.Error(); err != nil {
		return fmt.Errorf("csv flush: %w", err)
	}

	verb := "appended to"
	if fresh {
		verb = "written to"
	}
	fmt.Printf("\ncsv: %d rows %s %s\n", len(c.rows), verb, c.path)
	return nil
}

// csvIsFresh reports whether the file needs a header, and refuses one whose header is a different
// set of columns.
func csvIsFresh(path string) (bool, error) {
	f, err := os.Open(path)
	if os.IsNotExist(err) {
		return true, nil
	}
	if err != nil {
		return false, fmt.Errorf("csv: %w", err)
	}
	defer f.Close()

	head, err := csv.NewReader(f).Read()
	if err == io.EOF {
		return true, nil // the file exists but is empty
	}
	if err != nil {
		return false, fmt.Errorf("csv: reading the header of %s: %w", path, err)
	}
	if !slices.Equal(head, csvHeader) {
		return false, fmt.Errorf("csv: %s has %d columns starting %q, this run writes %d starting %q; append would corrupt it",
			path, len(head), strings.Join(head[:min(3, len(head))], ","), len(csvHeader), strings.Join(csvHeader[:3], ","))
	}
	return false, nil
}

func itoa(v int) string { return strconv.Itoa(v) }

func ftoa(v float64, prec int) string {
	if math.IsInf(v, 0) || math.IsNaN(v) {
		return ""
	}
	return strconv.FormatFloat(v, 'f', prec, 64)
}

// negLog2 turns a distance into a precision in bits. An exact hit has infinite precision and is
// left empty rather than written as +Inf, which no reader parses the same way.
func negLog2(d float64) float64 {
	if d <= 0 {
		return math.Inf(1)
	}
	return -math.Log2(d)
}
