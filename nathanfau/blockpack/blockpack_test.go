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

var (
	logNFlag = flag.Int("logn", 12, "CI ring log-degree (MaxSlots = 2^logn, capacity = 2^(logn-1) AES blocks)")
	// 13 is the smallest that fits one state; -whlogn 14 or more gives a capacity above 1, which
	// is what exercises the interleaving stride.
	whLogNFlag = flag.Int("whlogn", 13, "standard ring log-degree of the WH packing (capacity = 2^(logn-1)/4096 AES blocks)")
)

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

// TestWHPackRoundtrip is TestPackRoundtrip for the Walsh-Hadamard packing: the standard ring this
// time, since the WH circuit lives there, and one ciphertext instead of 64.
func TestWHPackRoundtrip(t *testing.T) {
	logN := *whLogNFlag

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            logN,
		LogQ:            []int{45},
		LogDefaultScale: 40,
	})
	if err != nil {
		t.Fatalf("NewParametersFromLiteral (logN=%d): %v", logN, err)
	}

	slots := params.MaxSlots()
	if slots != 1<<(logN-1) {
		t.Fatalf("MaxSlots = %d, want %d (the standard ring has N/2 complex slots)", slots, 1<<(logN-1))
	}
	capacity := WHCapacity(slots)
	if capacity == 0 {
		t.Fatalf("logN=%d gives %d slots, one WH state needs %d: use -whlogn 13 or more", logN, slots, WHSlotsPerBlock)
	}
	t.Logf("logN=%d  MaxSlots=%d  capacity=%d blocks  ciphertexts=1  stride=%d", logN, slots, capacity, capacity)

	sk := rlwe.NewKeyGenerator(params).GenSecretKeyNew()
	ecd := ckks.NewEncoder(params)
	enc := rlwe.NewEncryptor(params, sk)
	dec := rlwe.NewDecryptor(params, sk)

	seed := time.Now().UnixNano()
	t.Logf("random seed = %d", seed)
	rng := rand.New(rand.NewSource(seed))

	blocks := WHRandomBlocks(slots, rng)

	ct, err := WHEncrypt(params, ecd, enc, blocks, params.MaxLevel())
	if err != nil {
		t.Fatalf("WHEncrypt: %v", err)
	}

	if err := WHCheck(params, ecd, dec, ct, blocks); err != nil {
		t.Fatalf("WHEncrypt->WHDecrypt (seed=%d): %v", seed, err)
	}
	t.Logf("OK: %d/%d blocks match after WHEncrypt->WHDecrypt", len(blocks), len(blocks))
}
