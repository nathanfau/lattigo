package transciphering

import (
	"encoding/binary"
	"encoding/csv"
	"flag"
	"fmt"
	"hash/fnv"
	"io"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"slices"
	"strconv"
	"strings"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
)

// One row per operation of a TestAES run, APPENDED to -csv.
//
// Every row repeats the run's parameters, and carries the run's start time. That is deliberate:
// successive runs can share one file, a comparison is a group-by on run_ts or on whichever
// parameter changed, and no row is ever ambiguous about which run produced it.
//
//	go test ./nathanfau/transciphering/ -run '^TestAES$' -v -timeout 0 -seed 42 -csv runs/aes.csv
//	go test ./nathanfau/transciphering/ -run '^TestAES$' -v -timeout 0 -seed 42 -cleanextract -csv runs/aes.csv

var csvFlag = flag.String("csv", "", "append to this file: parameters, timing, and the precision pooled over all 64 ciphertexts. TestAES and TestTransciphering write one row per operation (same columns), TestAESBench a single row per run; the two carry different columns, so give them different files")

// csvRunHeader is what every row repeats: the run, and everything that defines what it ran on.
// csvRowHeader is what each operation adds. Splitting them lets newAESCSV check the widths match
// BEFORE the run starts, rather than discovering it when the rows are written at the end.
var csvRunHeader = []string{
	// which run, from which tree, on which machine
	"run_ts", "host", "gomaxprocs", "go_version", "git_commit",
	// what was run
	"logn", "k", "slots", "blocks", "seed", "rounds",
	"subbytes", "xor", "clean", "cleandepth", "clean_scale", "refresh_scale", "place", "extract", "extractlv",
	// the moduli chain
	"logscale", "primes", "logq0", "logq", "logp", "logqp",
	"dnum", "digit_bits", "marge", "q_sizes", "p_sizes", "chain_id",
	// the secret and the bootstrapping circuit
	"h", "h_tilde", "s2c_levels", "c2s_levels", "mod1_type", "mod1_deg", "mod1_k", "logmsgratio",
	// what the keys cost
	"galois_keys", "key_gb", "keygen_ms", "dft_ms",
	// where the pipeline sits on the chain
	"refresh_lv", "ark_lv",
}

var csvRowHeader = []string{
	// where in the run
	"seq", "round", "step",
	// what it cost
	"ms", "elapsed_ms",
	// Where the REFRESH spent its ms, stage by stage and in pipeline order; empty on every other
	// step. Every stage that is timed has its own column; `ms_reste` is what none of them claimed --
	// `ms` minus their sum, so the twelve ALWAYS reconcile and an unmeasured cost shows up there
	// instead of hiding inside a bucket.
	"ms_bitpack", "ms_conv", "ms_stc", "ms_scaledown", "ms_modup", "ms_cts", "ms_combine",
	"ms_evalmod", "ms_extract", "ms_sr", "ms_recombine", "ms_reste",
	// what came out
	"level", "prec_avg", "prec_min", "worst_err", "slots_pooled", "bit_err", "blocks_wrong",
}

var csvHeader = append(append([]string{}, csvRunHeader...), csvRowHeader...)

// aesCSV accumulates the rows and writes them once, at the end: a run that dies mid-way still
// leaves nothing half-written to misread later.
type aesCSV struct {
	path    string
	header  []string // the columns this file carries, run + row
	run     []string // the parameter columns, identical on every row
	rows    [][]string
	seq     int
	elapsed time.Duration // the operations so far, which is NOT the wall clock: see add
}

// newAESCSV builds the per-operation recorder of TestAES.
func newAESCSV(path string, ctx *Context, cfg Config, logN, k, blocks, rounds int, seed int64) (*aesCSV, error) {
	if path == "" {
		return nil, nil
	}
	if err := csvPrepare(path, csvHeader); err != nil {
		return nil, err
	}
	run, err := csvRunValues(ctx, cfg, logN, k, blocks, rounds, seed)
	if err != nil {
		return nil, err
	}
	return &aesCSV{path: path, header: csvHeader, run: run}, nil
}

// csvPrepare checks, before anything is computed, that the file can be appended to: finding out at
// the END of a ten-minute run that its header belongs to an older set of columns would be a poor
// way to learn it. An empty path is "no csv" and passes.
func csvPrepare(path string, header []string) error {
	if path == "" {
		return nil
	}
	if dir := filepath.Dir(path); dir != "." {
		if err := os.MkdirAll(dir, 0o755); err != nil {
			return fmt.Errorf("csv: %w", err)
		}
	}
	_, err := csvIsFresh(path, header)
	return err
}

// csvRunValues is what every row of a run repeats: the run itself, and everything that defines what
// it ran on. It needs the context, so it can only be built once the keys are generated.
func csvRunValues(ctx *Context, cfg Config, logN, k, blocks, rounds int, seed int64) ([]string, error) {
	p, btp := ctx.Params, ctx.Eval.Parameters
	Q, P := p.Q(), p.P()
	mod1P := btp.Mod1ParametersLiteral
	digit := widestDigit(Q, len(P))
	run := []string{
		time.Now().Format(time.RFC3339), hostname(), itoa(runtime.GOMAXPROCS(0)), runtime.Version(), gitCommit(),

		itoa(logN), itoa(k), itoa(ctx.Sw.CiP.MaxSlots()), itoa(blocks), strconv.FormatInt(seed, 10), itoa(rounds),

		sboxName(), cfg.Xor.String(), cfg.Clean.String(), itoa(cfg.Clean.Depth()), cleanScaleName(ctx), refreshScaleName(ctx),
		cfg.Place.String(), extractName(cfg), itoa(cfg.ExtractLevels(k)),

		itoa(p.LogDefaultScale()), itoa(len(Q)), itoa(sizeOf(Q[0])),
		ftoa(p.LogQ(), 1), ftoa(p.LogP(), 1), ftoa(p.LogQP(), 1),

		itoa(p.BaseRNSDecompositionVectorSize(p.MaxLevel(), p.MaxLevelP())),
		ftoa(digit, 1), ftoa(p.LogP()-digit, 1), chainSizes(Q), chainSizes(P), chainID(Q, P),

		itoa(p.XsHammingWeight()), itoa(btp.EphemeralSecretWeight),
		joinInts(btp.SlotsToCoeffsParameters.Levels), joinInts(btp.CoeffsToSlotsParameters.Levels),
		mod1Name(mod1P.Mod1Type), itoa(mod1P.Mod1Degree), itoa(mod1P.K), itoa(mod1P.LogMessageRatio),

		itoa(len(btp.GaloisElements(p))),
		ftoa(float64(ctx.BtpKeyBytes+ctx.SwKeyBytes)/(1<<30), 3),
		ftoa(millis(ctx.BtpKeyGen+ctx.SwKeyGen), 1), ftoa(millis(ctx.EvalSetup), 1),

		itoa(ctx.RefreshLv), itoa(ctx.ARKLv),
	}
	if len(run) != len(csvRunHeader) {
		return nil, fmt.Errorf("csv: %d run values for %d run columns", len(run), len(csvRunHeader))
	}
	return run, nil
}

// add records one operation. prec.WorstBit is a distance, so prec_min is its -log2: the MINIMUM
// precision over the pool, the number that decides whether a bit is about to flip.
//
// elapsed_ms is the running sum of the operations, i.e. what the PIPELINE has cost when this
// precision is reached -- not the wall clock, which also carries the oracle traces between the
// operations. That is the x of a precision-against-cost plot: two configurations that spend a
// different number of operations are only comparable on the time they spend.
// add records one operation. tm is the stage breakdown of a Refresh; pass the zero value for every
// other step and the seven columns come out empty rather than zero -- an empty cell reads as "this
// step has no stages", a zero would read as "it spent no time there".
func (c *aesCSV) add(round int, step string, dur time.Duration, prec precSummary, tm RefreshTimings) {
	if c == nil {
		return
	}
	c.seq++
	c.elapsed += dur
	row := append([]string{}, c.run...)
	row = append(row,
		itoa(c.seq), itoa(round), step,
		ftoa(float64(dur.Microseconds())/1000, 3),
		ftoa(float64(c.elapsed.Microseconds())/1000, 3))
	row = append(row, stageCols(tm, dur)...)
	row = append(row,
		itoa(prec.Level), ftoa(prec.AvgPrec, 3), ftoa(negLog2(prec.WorstBit), 3), ftoa(prec.WorstBit, 6),
		itoa(prec.Slots), ftoa(prec.BitErr, 6), itoa(prec.Wrong))
	c.rows = append(c.rows, row)
}

// stageCols renders the stage columns, or as many empty cells when nothing was measured. The last
// one is the remainder: dur minus everything the stages claimed, so the columns always add back up
// to `ms`. A stage nobody times lands there and is visible, instead of inflating a neighbour.
func stageCols(tm RefreshTimings, dur time.Duration) []string {
	if tm == (RefreshTimings{}) {
		return make([]string, 12)
	}
	stages := []time.Duration{
		tm.BitPack, tm.Conv, tm.STC, tm.ScaleDown, tm.ModUp, tm.CTS, tm.Combine,
		tm.EvalMod, tm.Extract, tm.SR, tm.Recombine,
	}
	out := make([]string, 0, len(stages)+1)
	rest := dur
	for _, d := range stages {
		out = append(out, millisOf(d))
		rest -= d
	}
	return append(out, millisOf(rest))
}

func millisOf(d time.Duration) string { return ftoa(float64(d.Microseconds())/1000, 3) }

// write APPENDS to the file, so successive runs accumulate: point every run at the same -csv and
// the comparison is a group-by on run_ts, or on whichever parameter column changed. The header goes
// in only when the file is new, and an existing file whose header is not this one is refused rather
// than corrupted -- that happens when the columns here change under an older file.
func (c *aesCSV) write() error {
	if c == nil {
		return nil
	}

	for i, r := range c.rows {
		if len(r) != len(c.header) {
			return fmt.Errorf("csv: row %d has %d fields for %d columns", i, len(r), len(c.header))
		}
	}

	fresh, err := csvIsFresh(c.path, c.header)
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
		if err := w.Write(c.header); err != nil {
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
func csvIsFresh(path string, header []string) (bool, error) {
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
	if !slices.Equal(head, header) {
		return false, fmt.Errorf("csv: %s has %d columns starting %q, this run writes %d starting %q; append would corrupt it",
			path, len(head), strings.Join(head[:min(3, len(head))], ","), len(header), strings.Join(header[:3], ","))
	}
	return false, nil
}

// cleanScaleName says where the round cleaning lands: "input" (its input's scale, the historical
// behaviour) or "fixed" (the default scale, -cleanfixed or -sbexact).
func cleanScaleName(ctx *Context) string {
	if ctx.CleanFixedScale || ctx.SBoxExact {
		return "fixed"
	}
	return "input"
}

// refreshScaleName says how the refresh hands its bits back on Canon: "relabel" (the label is
// overwritten, the historical behaviour) or "exact" (the extraction lands there, -refreshcanon).
func refreshScaleName(ctx *Context) string {
	if ctx.RefreshCanon {
		return "exact"
	}
	return "relabel"
}

func itoa(v int) string { return strconv.Itoa(v) }

func millis(d time.Duration) float64 { return float64(d.Microseconds()) / 1000 }

func joinInts(v []int) string {
	s := make([]string, len(v))
	for i, x := range v {
		s[i] = itoa(x)
	}
	return strings.Join(s, " ")
}

// sizeOf is a prime's size in bits, ROUNDED: an NTT-friendly prime sits a hair above its power of
// two, so its bit length is one more than the size it was asked for and reads as another prime.
func sizeOf(q uint64) int { return int(math.Round(math.Log2(float64(q)))) }

func chainSizes(mods []uint64) string {
	s := make([]string, len(mods))
	for i, q := range mods {
		s[i] = itoa(sizeOf(q))
	}
	return strings.Join(s, " ")
}

// widestDigit is the largest key-switching digit, len(P) consecutive Q primes, measured on the
// primes themselves: summing the rounded sizes under-reports it by a few bits, and those are the
// bits of margin against the key-switch noise.
func widestDigit(Q []uint64, nbPi int) (best float64) {
	for i := 0; i+nbPi <= len(Q); i++ {
		var s float64
		for j := i; j < i+nbPi; j++ {
			s += math.Log2(float64(Q[j]))
		}
		best = math.Max(best, s)
	}
	return best
}

// chainID identifies the primes themselves, which their sizes do not: change the size of any one
// of them and the generator hands a DIFFERENT prime to every level below, which alone moves the
// precision by bits. Two rows share a chain_id iff they ran on the same chain.
func chainID(Q, P []uint64) string {
	h := fnv.New64a()
	var b [8]byte
	for _, m := range append(append([]uint64{}, Q...), P...) {
		binary.LittleEndian.PutUint64(b[:], m)
		h.Write(b[:])
	}
	return fmt.Sprintf("%08x", uint32(h.Sum64()))
}

func mod1Name(t mod1.Type) string {
	switch t {
	case mod1.CosDiscrete:
		return "CosDiscrete"
	case mod1.SinContinuous:
		return "SinContinuous"
	case mod1.CosContinuous:
		return "CosContinuous"
	}
	return itoa(int(t))
}

func hostname() string {
	h, err := os.Hostname()
	if err != nil {
		return ""
	}
	return h
}

// gitCommit is the tree the run came from, with a trailing '+' when it was dirty, and empty when
// git cannot say. A file that accumulates runs over weeks is unreadable without it.
func gitCommit() string {
	out, err := exec.Command("git", "rev-parse", "--short", "HEAD").Output()
	if err != nil {
		return ""
	}
	rev := strings.TrimSpace(string(out))
	if st, err := exec.Command("git", "status", "--porcelain").Output(); err == nil && len(strings.TrimSpace(string(st))) > 0 {
		rev += "+"
	}
	return rev
}

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
