package algo1

import "time"

// Timings accumulates where Algorithm 1 spent its time, stage by stage, so a run can say what the
// refresh actually costs rather than only how long it took in total. Durations ADD UP across calls:
// the refresh runs ExtractTo sixteen times and Resume once, and one Timings collects the lot.
//
// The plain ExtractTo and Resume keep their signatures and time nothing; the Timed variants take a
// *Timings. A nil one is legal and turns the accounting off, which is what the callers that only
// want the result pass.
type Timings struct {
	STC       time.Duration // SlotsToCoeffs
	ScaleDown time.Duration // the rescaling that precedes the mod raise
	ModUp     time.Duration // the mod raise itself
	CTS       time.Duration // CoeffsToSlots
	Combine   time.Duration // the CombineReIm that closes step1
	Conv      time.Duration // the CI <-> Std conversions Algorithm 1 does itself
	EvalMod   time.Duration // EvalCos and EvalSin, plus Resume's double-angle squarings
	Extract   time.Duration // extractExp
}

// since adds the time elapsed since t0 to d. It exists so the instrumented call sites read as one
// line each and never branch: extractTo and resume normalise a nil *Timings to a throwaway one, so
// d is always a real field.
func since(d *time.Duration, t0 time.Time) { *d += time.Since(t0) }
