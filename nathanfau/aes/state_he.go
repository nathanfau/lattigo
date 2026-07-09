package aes

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// EncStateRepl encrypts a cleartext state bit-sliced, each bit replicated over all slots.
func EncStateRepl(params ckks.Parameters, ecd *ckks.Encoder, enc *rlwe.Encryptor, state [16]byte, level int) (st StateHE, err error) {
	slots := params.MaxSlots()
	for b := 0; b < 16; b++ {
		for i := 0; i < 8; i++ {
			bit := float64((state[b] >> uint(i)) & 1)
			vals := make([]float64, slots)
			for s := range vals {
				vals[s] = bit // replicate over all slots
			}
			pt := ckks.NewPlaintext(params, level)
			if err = ecd.Encode(vals, pt); err != nil {
				return st, fmt.Errorf("EncStateRepl byte %d bit %d: %w", b, i, err)
			}
			ct := ckks.NewCiphertext(params, 1, level)
			if err = enc.Encrypt(pt, ct); err != nil {
				return st, fmt.Errorf("EncStateRepl encrypt byte %d bit %d: %w", b, i, err)
			}
			st[b][i] = ct
		}
	}
	return st, nil
}

// DecStateBatch decrypts a bit-sliced state into one cleartext state per slot (bit threshold 0.5).
func DecStateBatch(params ckks.Parameters, ecd *ckks.Encoder, dec *rlwe.Decryptor, st StateHE) ([][16]byte, error) {
	slots := params.MaxSlots()
	out := make([][16]byte, slots)
	buf := make([]complex128, slots)
	for b := 0; b < 16; b++ {
		for i := 0; i < 8; i++ {
			if err := ecd.Decode(dec.DecryptNew(st[b][i]), buf); err != nil {
				return nil, fmt.Errorf("DecStateBatch byte %d bit %d: %w", b, i, err)
			}
			for s := 0; s < slots; s++ {
				if real(buf[s]) >= 0.5 { // bit in {0,1}, threshold 0.5
					out[s][b] |= 1 << uint(i)
				}
			}
		}
	}
	return out, nil
}
