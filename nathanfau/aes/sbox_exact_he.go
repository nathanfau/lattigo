package aes

// The S-box with EXACT scale alignment, side by side with sbox_he.go: same ANF, same balanced build,
// same laziness predicates, same levels.
//
// Why. sboxSum adds monomials built at different depths. Each rescale divides by a prime
// q = Delta(1+eta), so the monomials reach the summation level on scales that differ by products of
// (1+eta); FlattenLevels only drops levels, and lattigo's Add realigns two scales by their ratio
// TRUNCATED to an integer, i.e. adds them as if that ratio were 1. Weighted by the integer ANF
// (sum|c| ~ 2^8.7), that is the floor SubBytes sits on: about 2^-15 at logN 12 and 2^-9 at logN 16,
// where the primes, all = 1 mod 4N, lie further from their power of two.
//
// How. Every ciphertext carries, at each level l, one exact scale T_l: T_top is the scale of the
// input bits and T_{l-1} = T_l^2 / q_l, the scale a product of two T_l ciphertexts lands on after its
// rescale. A ciphertext that has to go down to a lower level consumes the LAST level of its drop
// with a multiplication by the integer round(r . q), r the exact ratio of the target scale to its
// own, followed by a rescale by q: its value is unchanged and it lands on the target scale. Those
// are levels DropLevel would have thrown away, so the circuit keeps the original's levels: no prime
// and no key switch more, one rescale more per drop. Every label is then pinned to its exact value,
// so each Add compares equal scales and never realigns.

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// SubByteExact applies the exactly aligned S-box variant named by version (1..3), the counterpart
// of SubByte.
func (a *Evaluator) SubByteExact(inp ByteHE, version int) (ByteHE, error) {
	switch version {
	case 1:
		return a.SubByteV1Exact(inp)
	case 2:
		return a.SubByteV2Exact(inp)
	case 3:
		return a.SubByteV3Exact(inp)
	default:
		return ByteHE{}, fmt.Errorf("SubByteExact: unknown version %d (want 1..3)", version)
	}
}

// SubByteV1Exact, SubByteV2Exact and SubByteV3Exact are SubByteV1, V2 and V3 with exact alignment:
// the same relin and main rescale counts, plus one alignment rescale per level drop.
func (a *Evaluator) SubByteV1Exact(inp ByteHE) (ByteHE, error) { return a.subByteExact(inp, 0) }
func (a *Evaluator) SubByteV2Exact(inp ByteHE) (ByteHE, error) { return a.subByteExact(inp, 4) }
func (a *Evaluator) SubByteV3Exact(inp ByteHE) (ByteHE, error) { return a.subByteExact(inp, 2) }

func (a *Evaluator) subByteExact(inp ByteHE, lazyMinDeg int) (ByteHE, error) {
	l, err := newLadder(a.eval, inp)
	if err != nil {
		return ByteHE{}, err
	}
	s := anf()
	isLazy := s.lazyAbove(lazyMinDeg)
	mono, err := l.buildMonomials(inp, s, isLazy)
	if err != nil {
		return ByteHE{}, err
	}
	return l.sboxSum(mono, s)
}

// ladder is the exact scale of every level at and below the S-box input.
type ladder struct {
	eval   *ckks.Evaluator
	params ckks.Parameters
	scale  map[int]rlwe.Scale // T_l, the scale of a degree-1 ciphertext at level l

	// aligned keeps every ciphertext already brought down, per source and level: a factor that
	// several products use at the same level, and that the sum uses again, is aligned once.
	aligned map[alignKey]*rlwe.Ciphertext
}

type alignKey struct {
	src   *rlwe.Ciphertext
	level int
}

// newLadder builds the ladder from the input bits, which must share one level and one scale. As in
// pin, a scale within 2^-40 of bit 0's is the same scale up to the rounding of lattigo's label
// arithmetic (the polynomial evaluator of the cleaning leaves such gaps): it is pinned to bit 0's.
func newLadder(eval *ckks.Evaluator, inp ByteHE) (*ladder, error) {
	params := *eval.GetParameters()
	if n := params.LevelsConsumedPerRescaling(); n != 1 {
		return nil, fmt.Errorf("exact S-box: %d primes per rescale, only 1 is handled", n)
	}
	top, T := inp[0].Level(), inp[0].Scale
	for j, ct := range inp {
		if ct.Level() != top || ct.Degree() != 1 || !ct.Scale.InDelta(T, 40) {
			return nil, fmt.Errorf("exact S-box: input bit %d at level %d, degree %d, scale 2^%.6f (off by 2^-%.1f); bit 0 at level %d, scale 2^%.6f: the bits must share one level and one scale",
				j, ct.Level(), ct.Degree(), ct.Scale.Log2(), ct.Scale.Log2Delta(T), top, T.Log2())
		}
		ct.Scale = T
	}
	l := &ladder{eval: eval, params: params, scale: map[int]rlwe.Scale{top: T}, aligned: map[alignKey]*rlwe.Ciphertext{}}
	for lv := top; lv > 0; lv-- {
		l.scale[lv-1] = l.scale[lv].Mul(l.scale[lv]).Div(utils.RescalingFactor(params, lv))
	}
	return l, nil
}

// want is the exact scale of a ciphertext of this degree at this level: T_l, or T_l^2 for a product
// not relinearized yet.
func (l *ladder) want(level, degree int) (rlwe.Scale, error) {
	t, ok := l.scale[level]
	if !ok {
		return rlwe.Scale{}, fmt.Errorf("exact S-box: level %d is above the input", level)
	}
	switch degree {
	case 1:
		return t, nil
	case 2:
		return t.Mul(t), nil
	}
	return rlwe.Scale{}, fmt.Errorf("exact S-box: degree %d ciphertext", degree)
}

// pin checks that ct carries its exact scale up to the rounding of lattigo's own label arithmetic,
// then writes that exact value, so every later Add compares equal labels.
func (l *ladder) pin(ct *rlwe.Ciphertext) error {
	want, err := l.want(ct.Level(), ct.Degree())
	if err != nil {
		return err
	}
	if !ct.Scale.InDelta(want, 40) {
		return fmt.Errorf("exact S-box: level %d degree %d, scale 2^%.6f is off its exact value by 2^-%.1f",
			ct.Level(), ct.Degree(), ct.Scale.Log2(), ct.Scale.Log2Delta(want))
	}
	ct.Scale = want
	return nil
}

// dropTo returns a copy of ct at the given level, same value, exact scale. The levels above the
// last one are dropped; the last one is spent on round(r . q) x ct and a rescale by q. The result is
// cached, so the same ciphertext is aligned to the same level only once; callers get their own copy.
func (l *ladder) dropTo(ct *rlwe.Ciphertext, level int) (*rlwe.Ciphertext, error) {
	key := alignKey{ct, level}
	if c, ok := l.aligned[key]; ok {
		return c.CopyNew(), nil
	}
	out, err := l.align(ct, level)
	if err != nil {
		return nil, err
	}
	l.aligned[key] = out
	return out.CopyNew(), nil
}

// align is dropTo without the cache.
func (l *ladder) align(ct *rlwe.Ciphertext, level int) (*rlwe.Ciphertext, error) {
	out := ct.CopyNew()
	d := out.Level() - level
	switch {
	case d < 0:
		return nil, fmt.Errorf("exact S-box: cannot raise level %d to %d", out.Level(), level)
	case d == 0:
		return out, l.pin(out)
	case d > 1:
		l.eval.DropLevel(out, d-1)
	}
	want, err := l.want(level, out.Degree())
	if err != nil {
		return nil, err
	}
	if want.Equal(out.Scale) {
		l.eval.DropLevel(out, 1)
		return out, nil
	}

	// K = round(want / scale . q): an INTEGER constant, which lattigo multiplies without touching
	// the label, so the rescale that follows divides the polynomial by q exactly as it divides the
	// label, and the polynomial ends up carrying value x want.
	q := utils.RescalingFactor(l.params, level+1)
	k := new(big.Float).SetPrec(256).Quo(&want.Value, &out.Scale.Value)
	k.Mul(k, &q.Value)
	k.Add(k, big.NewFloat(0.5))
	K, _ := k.Int(nil)

	if err := l.eval.Mul(out, K, out); err != nil {
		return nil, fmt.Errorf("exact S-box: align Mul: %w", err)
	}
	if err := l.eval.Rescale(out, out); err != nil {
		return nil, fmt.Errorf("exact S-box: align Rescale: %w", err)
	}
	utils.Ops.Align++
	out.Scale = want
	return out, nil
}

// mul is utils.MulLeveled (or MulLeveledLazy) with exact drops: both operands meet at the lower
// level on its exact scale.
func (l *ladder) mul(a, b *rlwe.Ciphertext, lazy bool) (*rlwe.Ciphertext, error) {
	lv := min(a.Level(), b.Level())
	p, err := l.dropTo(a, lv)
	if err != nil {
		return nil, err
	}
	q, err := l.dropTo(b, lv)
	if err != nil {
		return nil, err
	}
	if lazy {
		out, err := l.eval.MulNew(p, q)
		if err != nil {
			return nil, fmt.Errorf("exact S-box: lazy Mul: %w", err)
		}
		return out, l.pin(out)
	}
	if err = l.eval.MulRelin(p, q, p); err != nil {
		return nil, fmt.Errorf("exact S-box: MulRelin: %w", err)
	}
	utils.Ops.Relin++
	if err = l.eval.Rescale(p, p); err != nil {
		return nil, fmt.Errorf("exact S-box: Rescale: %w", err)
	}
	utils.Ops.Rescale++
	return p, l.pin(p)
}

// buildMonomials is Evaluator.buildMonomials on exact multiplications.
func (l *ladder) buildMonomials(inp ByteHE, s *sboxANF, isLazy lazyFn) (map[int]*rlwe.Ciphertext, error) {
	mono := make(map[int]*rlwe.Ciphertext, len(s.needed)+8)
	for j := 0; j < 8; j++ {
		mono[1<<uint(j)] = inp[j]
	}
	var build func(S int) (*rlwe.Ciphertext, error)
	build = func(S int) (*rlwe.Ciphertext, error) {
		if m, ok := mono[S]; ok {
			return m, nil
		}
		lo, hi := balancedSplit(S)
		lc, err := build(lo)
		if err != nil {
			return nil, err
		}
		hc, err := build(hi)
		if err != nil {
			return nil, err
		}
		m, err := l.mul(lc, hc, isLazy(S))
		if err != nil {
			return nil, fmt.Errorf("buildMonomials S=%d: %w", S, err)
		}
		mono[S] = m
		return m, nil
	}
	for _, S := range s.needed {
		if _, err := build(S); err != nil {
			return nil, err
		}
	}
	return mono, nil
}

// sboxSum is Evaluator.sboxSum on exact scales: every monomial goes down to the summation level on
// that level's exact scale (T_m, or T_m^2 for a lazy one) instead of being dropped, so the eager and
// the lazy accumulators each add equal labels, and the eager one goes down to the lazy one's level
// the same way before the two are combined.
func (l *ladder) sboxSum(mono map[int]*rlwe.Ciphertext, s *sboxANF) (out ByteHE, err error) {
	m := -1
	for _, S := range s.needed {
		if lv := mono[S].Level(); m < 0 || lv < m {
			m = lv
		}
	}
	terms := make(map[int]*rlwe.Ciphertext, len(s.needed))
	for _, S := range s.needed {
		if terms[S], err = l.dropTo(mono[S], m); err != nil {
			return out, fmt.Errorf("sboxSum S=%d: %w", S, err)
		}
	}
	eval := l.eval

	addTerm := func(accPtr **rlwe.Ciphertext, S, c int) error {
		t := terms[S]
		acc := *accPtr
		switch {
		case acc == nil:
			acc = t.CopyNew()
			if c != 1 {
				if e := eval.Mul(acc, c, acc); e != nil {
					return e
				}
			}
		case c == 1:
			if e := eval.Add(acc, t, acc); e != nil {
				return e
			}
		case c == -1:
			if e := eval.Sub(acc, t, acc); e != nil {
				return e
			}
		default:
			tc := t.CopyNew()
			if e := eval.Mul(tc, c, tc); e != nil {
				return e
			}
			if e := eval.Add(acc, tc, acc); e != nil {
				return e
			}
		}
		*accPtr = acc
		return nil
	}

	lowest := m
	for i := 0; i < 8; i++ {
		var accE, accL *rlwe.Ciphertext
		for _, S := range s.needed {
			c := s.coeffs[i][S]
			if c == 0 {
				continue
			}
			acc, kind := &accE, "eager"
			if terms[S].Degree() == 2 {
				acc, kind = &accL, "lazy"
			}
			if err = addTerm(acc, S, c); err != nil {
				return out, fmt.Errorf("sboxSum bit %d %s S=%d: %w", i, kind, S, err)
			}
		}
		if accL != nil {
			if err = eval.Relinearize(accL, accL); err != nil {
				return out, fmt.Errorf("sboxSum bit %d relin: %w", i, err)
			}
			utils.Ops.Relin++
			if err = eval.Rescale(accL, accL); err != nil {
				return out, fmt.Errorf("sboxSum bit %d rescale: %w", i, err)
			}
			utils.Ops.Rescale++
			if err = l.pin(accL); err != nil {
				return out, fmt.Errorf("sboxSum bit %d: %w", i, err)
			}
		}
		var acc *rlwe.Ciphertext
		switch {
		case accE != nil && accL != nil:
			if accE, err = l.dropTo(accE, accL.Level()); err != nil {
				return out, fmt.Errorf("sboxSum bit %d combine: %w", i, err)
			}
			if err = eval.Add(accE, accL, accE); err != nil {
				return out, fmt.Errorf("sboxSum bit %d combine: %w", i, err)
			}
			acc = accE
		case accE != nil:
			acc = accE
		case accL != nil:
			acc = accL
		default:
			return out, fmt.Errorf("sboxSum bit %d: empty sum", i)
		}
		if a0 := s.coeffs[i][0]; a0 != 0 {
			if err = eval.Add(acc, a0, acc); err != nil {
				return out, fmt.Errorf("sboxSum bit %d add const: %w", i, err)
			}
		}
		out[i] = acc
		lowest = min(lowest, acc.Level())
	}

	// The 8 output bits leave on one level and one scale, which is what the next circuit expects.
	for i := range out {
		if out[i], err = l.dropTo(out[i], lowest); err != nil {
			return out, fmt.Errorf("sboxSum bit %d out: %w", i, err)
		}
	}
	return out, nil
}
