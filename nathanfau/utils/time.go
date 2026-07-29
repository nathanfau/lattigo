package utils

// Timing helpers for instrumenting a circuit: they measure and report durations only, nothing
// CKKS-specific.

import (
	"fmt"
	"sort"
	"time"
)

// TimeSince prints the elapsed time since t0 under the given label + returns it
func TimeSince(label string, t0 time.Time) time.Duration {
	d := time.Since(t0)
	fmt.Printf("%-20s %s\n", label, d.Round(time.Microsecond))
	return d
}

// TimeStats prints min/max/mean/median (rounded to ms) over a list of per-step durations under the
// given label. Used to summarize a step's per-round timings across a multi-round run. If a positive
// reference duration is passed, it also prints the mean as a percentage of that reference
func TimeStats(name string, ds []time.Duration, ref ...time.Duration) {
	if len(ds) == 0 {
		fmt.Printf("  %-15s (no measure)\n", name)
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
	share := ""
	if len(ref) > 0 && ref[0] > 0 {
		share = fmt.Sprintf(" share=%5.1f%%", 100*float64(mean)/float64(ref[0]))
	}
	fmt.Printf("  %-15s min=%s max=%s mean=%s med=%s%s (n=%d)\n", name,
		cp[0].Round(time.Millisecond), cp[len(cp)-1].Round(time.Millisecond),
		mean.Round(time.Millisecond), med.Round(time.Millisecond), share, len(cp))
}
