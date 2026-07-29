package utils

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// WantToCplx converts a want vector ([]float64 or []complex128) to at most n complex128 values, so
// several can be concatenated for ckks.GetPrecisionStats. Returns nil for any other type.
func WantToCplx(want any, n int) []complex128 {
	switch w := want.(type) {
	case []float64:
		if len(w) < n {
			n = len(w)
		}
		out := make([]complex128, n)
		for i := 0; i < n; i++ {
			out[i] = complex(w[i], 0)
		}
		return out
	case []complex128:
		if len(w) < n {
			n = len(w)
		}
		return append([]complex128(nil), w[:n]...)
	default:
		return nil
	}
}

// BitStats measures how far the slots of one decoded plane sit from their reference values: the
// worst slot decides how much margin is left before a bit flips.
type BitStats struct {
	Worst     float64 // worst |value - reference| over the plane
	WorstSlot int
	Best      float64 // smallest such distance
	BestSlot  int
	AvgPrec   float64 // precision in bits averaged over the slots (mean of the -log2 distances)
	Wrong     int     // slots farther than Tol from their reference (0 when Tol <= 0)
	Slots     int     // slots measured
	Tol       float64 // tolerance applied, 0 when none was
}

// BitDistanceOf measures one ALREADY DECODED plane. want holds the reference values ([]float64 or
// []complex128, the distance being the real gap or the complex modulus accordingly).
// pass nil to measure each slot against the NEAREST bit instead, which says whether the values
// are clean bits without needing an oracle. A tol <= 0 skips the Wrong count.
func BitDistanceOf(got []complex128, want any, tol float64) BitStats {
	s := BitStats{Best: math.Inf(1), Tol: tol, WorstSlot: -1, BestSlot: -1}
	_, isCplx := want.([]complex128)
	ref := WantToCplx(want, len(got))
	if want != nil && ref == nil {
		return s
	}
	n := len(got)
	if ref != nil && len(ref) < n {
		n = len(ref)
	}
	s.Slots = n

	sumPrec, nPrec := 0.0, 0
	for i := 0; i < n; i++ {
		var d float64
		if ref == nil {
			// No oracle: distance to the nearest bit, 0 or 1.
			if v := real(got[i]); v >= 0.5 {
				d = math.Abs(v - 1)
			} else {
				d = math.Abs(v)
			}
		} else if diff := got[i] - ref[i]; isCplx {
			d = math.Hypot(real(diff), imag(diff))
		} else {
			d = math.Abs(real(diff))
		}
		if d > s.Worst {
			s.Worst, s.WorstSlot = d, i
		}
		if d < s.Best {
			s.Best, s.BestSlot = d, i
		}
		if tol > 0 && d > tol {
			s.Wrong++
		}
		if d > 0 {
			sumPrec += -math.Log2(d)
			nPrec++
		}
	}
	// An exactly-hit slot has infinite precision and is left out of the average.
	s.AvgPrec = math.Inf(1)
	if nPrec > 0 {
		s.AvgPrec = sumPrec / float64(nPrec)
	}
	return s
}

// BitDistanceCt decrypts and decodes ct, then measures it with BitDistanceOf.
func BitDistanceCt(ecd *ckks.Encoder, dec *rlwe.Decryptor, ct *rlwe.Ciphertext, want any, tol float64) (BitStats, error) {
	buf := make([]complex128, ct.Slots())
	if err := ecd.Decode(dec.DecryptNew(ct), buf); err != nil {
		return BitStats{}, fmt.Errorf("BitDistanceCt decode: %w", err)
	}
	return BitDistanceOf(buf, want, tol), nil
}

// BitDistance measures every ciphertext of cts with BitDistanceOf, one BitStats each. want[i] is
// the reference plane of cts[i]; a nil want measures against the nearest bit.
func BitDistance(ecd *ckks.Encoder, dec *rlwe.Decryptor, cts []*rlwe.Ciphertext, want [][]float64, tol float64) ([]BitStats, error) {
	if want != nil && len(want) != len(cts) {
		return nil, fmt.Errorf("BitDistance: %d ciphertexts for %d reference planes", len(cts), len(want))
	}
	out := make([]BitStats, len(cts))
	for i, ct := range cts {
		var ref any
		if want != nil {
			ref = want[i]
		}
		s, err := BitDistanceCt(ecd, dec, ct, ref, tol)
		if err != nil {
			return nil, fmt.Errorf("BitDistance ct %d: %w", i, err)
		}
		out[i] = s
	}
	return out, nil
}

// WorstOf folds per-ciphertext stats into one: the worst and best slots win, the counts add. It
// also returns the index of the ciphertext holding the worst slot (-1 if all is empty).
func WorstOf(all []BitStats) (BitStats, int) {
	agg := BitStats{Best: math.Inf(1), WorstSlot: -1, BestSlot: -1}
	worstCt := -1
	sumPrec, nPrec := 0.0, 0
	for i, s := range all {
		if s.Slots == 0 {
			continue
		}
		if s.Worst > agg.Worst {
			agg.Worst, agg.WorstSlot, worstCt = s.Worst, s.WorstSlot, i
		}
		if s.Best < agg.Best {
			agg.Best, agg.BestSlot = s.Best, s.BestSlot
		}
		if !math.IsInf(s.AvgPrec, 0) {
			sumPrec += s.AvgPrec * float64(s.Slots)
			nPrec += s.Slots
		}
		agg.Wrong += s.Wrong
		agg.Slots += s.Slots
		agg.Tol = s.Tol
	}
	agg.AvgPrec = math.Inf(1)
	if nPrec > 0 {
		agg.AvgPrec = sumPrec / float64(nPrec)
	}
	return agg, worstCt
}
