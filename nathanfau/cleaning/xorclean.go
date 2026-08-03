package cleaning

// The three ways to spend one cleaning and one XOR on a pair of bits. All three cost clean's depth
// plus one, so they are interchangeable in a level map; what changes is WHERE the cleaning lands.

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/aes"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// CleanFunc is the shape the three cleaning polynomials of this package share.
type CleanFunc func(ckks.Parameters, *ckks.Evaluator, *rlwe.Ciphertext) (*rlwe.Ciphertext, error)

// Kind names one of the three polynomials of cleaning.go.
type Kind int

const (
	Basic        Kind = iota // Cleaning, 2 levels, the pipeline default and the zero value
	Smoother                 // SmootherCleaning, 3 levels
	VerySmoother             // VerySmootherCleaning, 3 levels
)

func (k Kind) Depth() int {
	if k == Basic {
		return 2
	}
	return 3
}

func (k Kind) Func() CleanFunc {
	switch k {
	case Smoother:
		return SmootherCleaning
	case VerySmoother:
		return VerySmootherCleaning
	default:
		return Cleaning
	}
}

func (k Kind) String() string {
	switch k {
	case Smoother:
		return "SmootherCleaning (3 lv)"
	case VerySmoother:
		return "VerySmootherCleaning (3 lv)"
	default:
		return "Cleaning (2 lv)"
	}
}

// ParseKind reads the name a flag carries.
func ParseKind(s string) (Kind, error) {
	switch s {
	case "smoother", "s":
		return Smoother, nil
	case "cleaning", "basic", "c":
		return Basic, nil
	case "verysmoother", "vs":
		return VerySmoother, nil
	}
	return 0, fmt.Errorf("unknown cleaning %q, want \"smoother\", \"cleaning\" or \"verysmoother\"", s)
}

// RouteFunc is the shape the three routes below share.
type RouteFunc func(params ckks.Parameters, eval *ckks.Evaluator, clean CleanFunc, xor aes.XorFunc,
	x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error)

// Placement names WHERE the cleaning lands relative to the XOR, i.e. which of the three routes runs.
// The zero value is XorThenClean, the historical pipeline.
type Placement int

const (
	CleanAfter Placement = iota // XorThenClean, h(x XOR y)
	CleanBoth                   // CleanThenXor, h(x) XOR h(y)
	CleanOne                    // CleanOneThenXor, h(x) XOR y
)

func (p Placement) Func() RouteFunc {
	switch p {
	case CleanBoth:
		return CleanThenXor
	case CleanOne:
		return CleanOneThenXor
	default:
		return XorThenClean
	}
}

func (p Placement) String() string {
	switch p {
	case CleanBoth:
		return "CleanThenXor"
	case CleanOne:
		return "CleanOneThenXor"
	default:
		return "XorThenClean"
	}
}

// YDrop is how many levels BELOW x the second operand y has to sit for the route to line up.
func (p Placement) YDrop(cleanDepth int) int {
	if p == CleanOne {
		return cleanDepth
	}
	return 0
}

// ParsePlacement reads the name a flag carries.
func ParsePlacement(s string) (Placement, error) {
	switch s {
	case "after", "xorthenclean", "a":
		return CleanAfter, nil
	case "both", "cleanthenxor", "b":
		return CleanBoth, nil
	case "one", "cleanonethenxor", "o":
		return CleanOne, nil
	}
	return 0, fmt.Errorf("unknown placement %q, want \"after\", \"both\" or \"one\"", s)
}

// XorThenClean cleans the RESULT: h(x XOR y)
func XorThenClean(params ckks.Parameters, eval *ckks.Evaluator, clean CleanFunc, xor aes.XorFunc,
	x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	d, err := xor(x, y)
	if err != nil {
		return nil, fmt.Errorf("XorThenClean xor: %w", err)
	}
	out, err := clean(params, eval, d)
	if err != nil {
		return nil, fmt.Errorf("XorThenClean clean: %w", err)
	}
	return out, nil
}

// CleanThenXor cleans BOTH operands: h(x) XOR h(y)
func CleanThenXor(params ckks.Parameters, eval *ckks.Evaluator, clean CleanFunc, xor aes.XorFunc,
	x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	hx, err := clean(params, eval, x)
	if err != nil {
		return nil, fmt.Errorf("CleanThenXor h(x): %w", err)
	}
	hy, err := clean(params, eval, y)
	if err != nil {
		return nil, fmt.Errorf("CleanThenXor h(y): %w", err)
	}
	utils.AlignLevels(eval, hx, hy)

	out, err := xor(hx, hy)
	if err != nil {
		return nil, fmt.Errorf("CleanThenXor xor: %w", err)
	}
	return out, nil
}

// CleanOneThenXor cleans x only: h(x) XOR y
func CleanOneThenXor(params ckks.Parameters, eval *ckks.Evaluator, clean CleanFunc, xor aes.XorFunc,
	x, y *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	hx, err := clean(params, eval, x)
	if err != nil {
		return nil, fmt.Errorf("CleanOneThenXor h(x): %w", err)
	}
	yy := y.CopyNew()
	utils.AlignLevels(eval, hx, yy)

	out, err := xor(hx, yy)
	if err != nil {
		return nil, fmt.Errorf("CleanOneThenXor xor: %w", err)
	}
	return out, nil
}
