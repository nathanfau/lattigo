package bitbatching

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// bitPack packs the bit-ciphertexts bits[i] into a single ciphertext encrypting
// sum_i bits[i] * 2^i.
func BitPack(heEval *ckks.Evaluator, bits []*rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	ctPack := bits[0].CopyNew()
	tmp := bits[0].CopyNew()
	for i := 1; i < len(bits); i++ {
		if err := heEval.Mul(bits[i], math.Exp2(float64(i)), tmp); err != nil {
			return nil, fmt.Errorf("bitPack Mul i=%d: %w", i, err)
		}
		if err := heEval.Add(ctPack, tmp, ctPack); err != nil {
			return nil, fmt.Errorf("bitPack Add i=%d: %w", i, err)
		}
	}
	return ctPack, nil
}

// bitPackGroups splits bits into consecutive groups of k and packs each group with bitPack.
func BitPackGroups(heEval *ckks.Evaluator, bits []*rlwe.Ciphertext, k int) ([]*rlwe.Ciphertext, error) {
	if k < 1 {
		return nil, fmt.Errorf("bitPackGroups: k must be >= 1, got %d", k)
	}
	j := len(bits)
	if j%k != 0 {
		return nil, fmt.Errorf("bitPackGroups: %d ciphertexts not divisible by k=%d", j, k)
	}

	groups := make([]*rlwe.Ciphertext, j/k)
	for g := range groups {
		ctPack, err := BitPack(heEval, bits[g*k:(g+1)*k])
		if err != nil {
			return nil, fmt.Errorf("bitPackGroups group %d: %w", g, err)
		}
		groups[g] = ctPack
	}
	return groups, nil
}
