package blockpack

import (
	"flag"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var logNFlag = flag.Int("logn", 12, "CI ring log-degree (MaxSlots = 2^logn, capacity = 2^(logn-1) AES blocks)")

// TestPackRoundtrip checks the block-batched CI packing end to end at the chosen logN
func TestPackRoundtrip(t *testing.T) {
	logN := *logNFlag

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            logN,
		LogQ:            []int{45},
		LogDefaultScale: 40,
		RingType:        ring.ConjugateInvariant,
	})
	if err != nil {
		t.Fatalf("NewParametersFromLiteral (logN=%d): %v", logN, err)
	}

	slots := params.MaxSlots()
	if slots != 1<<logN {
		t.Fatalf("MaxSlots = %d, want %d (CI has N real slots)", slots, 1<<logN)
	}
	capacity := Capacity(params)
	t.Logf("logN=%d  MaxSlots=%d  capacity=%d blocks  ciphertexts=%d", logN, slots, capacity, 64)

	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	ecd := ckks.NewEncoder(params)
	enc := rlwe.NewEncryptor(params, sk)
	dec := rlwe.NewDecryptor(params, sk)

	seed := time.Now().UnixNano()
	t.Logf("random seed = %d", seed)
	rng := rand.New(rand.NewSource(seed))

	blocks := RandomBlocks(params, rng)

	p, err := Encrypt(params, ecd, enc, blocks, params.MaxLevel())
	if err != nil {
		t.Fatalf("Encrypt: %v", err)
	}

	if err := Check(params, ecd, dec, p, blocks); err != nil {
		t.Fatalf("Encrypt->Decrypt (seed=%d): %v", seed, err)
	}
	t.Logf("OK: %d/%d blocks match after Encrypt->Decrypt", len(blocks), len(blocks))
}
