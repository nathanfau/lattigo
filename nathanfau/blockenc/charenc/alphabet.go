package charenc

import "fmt"

type EncodingKind int

const (
	BRU EncodingKind = iota
	LBRU
	WH
	IND
)

func (k EncodingKind) String() string {
	switch k {
	case BRU:
		return "BRU"
	case LBRU:
		return "LBRU"
	case WH:
		return "WH"
	case IND:
		return "IND"
	default:
		return fmt.Sprintf("EncodingKind(%d)", int(k))
	}
}

type Alphabet struct {
	Modulus int
	Kind    EncodingKind
	// Omitted: only used by reduced IND.
	Omitted int
	// Generator: only used by LBRU (discrete-log basis of Z_t^*).
	Generator int
}

type BlockSpec struct {
	Alphabet Alphabet
	// Reduced selects the t-1 slot form vs the t slot form. LBRU interprets
	// this differently from the BRU/WH "omit the trivial character" rule.
	Reduced bool
	// Slots must match the encoding kind and Reduced flag; set by the codec.
	Slots int
}
