package utils

import (
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// RescalingFactor is what Rescale divides a ciphertext at this level by: the product of the primes
// it consumes, which is one of them unless the parameters rescale over several.
func RescalingFactor(params ckks.Parameters, level int) rlwe.Scale {
	ringQ := params.RingQ()
	s := rlwe.NewScale(1)
	for i := 0; i < params.LevelsConsumedPerRescaling(); i++ {
		s = s.Mul(rlwe.NewScale(ringQ.SubRings[level-i].Modulus))
	}
	return s
}
