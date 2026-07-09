package aes

import (
	"fmt"
	"testing"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// a fixed AES state used by every operation test.
var testState = [16]byte{
	0x19, 0x3d, 0xe3, 0xbe, 0xa0, 0xf4, 0xe2, 0x2b,
	0x9a, 0xc6, 0x8d, 0x2a, 0xe9, 0xf8, 0x48, 0x08,
}

func newAESContext(t *testing.T) (ckks.Parameters, *ckks.Encoder, *rlwe.Encryptor, *rlwe.Decryptor, *Evaluator) {
	t.Helper()
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            12,
		LogQ:            []int{55, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{61},
		LogDefaultScale: 45,
	})
	if err != nil {
		t.Fatalf("params: %v", err)
	}
	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	evk := rlwe.NewMemEvaluationKeySet(rlk)

	ecd := ckks.NewEncoder(params)
	enc := rlwe.NewEncryptor(params, pk)
	dec := rlwe.NewDecryptor(params, sk)
	ae := NewEvaluator(ckks.NewEvaluator(params, evk))

	debug.EncStd, debug.DecStd, debug.ParamsStd = ecd, dec, params

	return params, ecd, enc, dec, ae
}

func dbgState(label string, ae *Evaluator, params ckks.Parameters, clear [16]byte, st StateHE) {
	fmt.Printf(" cleartext %-4s = %x\n", label, clear)
	ct := st[0][0]
	debug.DbgSlotStd("ct_"+label+" bit[0][0] =", ct)
	debug.DbgChain("chain ct_"+label+" :", ae.eval, ct)
	want := make([]float64, params.MaxSlots())
	for s := range want {
		want[s] = float64(clear[0] & 1)
	}
	debug.PrintPrecStd("prec ct_"+label+" :", want, ct)
}

// checkState compares slot 0 of the decrypted state against the cleartext want.
func checkState(t *testing.T, params ckks.Parameters, ecd *ckks.Encoder, dec *rlwe.Decryptor, st StateHE, want [16]byte, op string) {
	t.Helper()
	got, err := DecStateBatch(params, ecd, dec, st)
	if err != nil {
		t.Fatalf("%s decrypt: %v", op, err)
	}
	if got[0] != want {
		t.Fatalf("%s mismatch:\n got  = %x\n want = %x", op, got[0], want)
	}
}

func TestAddRoundKey(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)
	rk := [16]byte{
		0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6,
		0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c,
	}

	st, err := EncStateRepl(params, ecd, enc, testState, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc state: %v", err)
	}
	rkHE, err := EncStateRepl(params, ecd, enc, rk, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc rk: %v", err)
	}

	want := testState
	AddRoundKey(want[:], rk[:])

	dbgState("in", ae, params, testState, st)
	out, err := ae.AddRoundKey(st, rkHE)
	if err != nil {
		t.Fatalf("AddRoundKey: %v", err)
	}
	dbgState("out", ae, params, want, out)

	checkState(t, params, ecd, dec, out, want, "AddRoundKey")
}

func TestSubBytesV1(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)

	st, err := EncStateRepl(params, ecd, enc, testState, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc state: %v", err)
	}

	want := testState
	SubBytes(want[:])

	dbgState("in", ae, params, testState, st)
	utils.ResetOps()
	var out StateHE
	for b := 0; b < 16; b++ {
		if out[b], err = ae.SubByteV1(st[b]); err != nil {
			t.Fatalf("SubByteV1 %d: %v", b, err)
		}
	}
	ops := utils.Ops
	t.Logf("SubBytesV1 (16 bytes): %d relin, %d rescale (%.1f relin, %.1f rescale per byte)",
		ops.Relin, ops.Rescale, float64(ops.Relin)/16, float64(ops.Rescale)/16)
	dbgState("out", ae, params, want, out)

	checkState(t, params, ecd, dec, out, want, "SubBytesV1")
}

func TestSubBytesV2(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)

	st, err := EncStateRepl(params, ecd, enc, testState, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc state: %v", err)
	}

	want := testState
	SubBytes(want[:])

	dbgState("in", ae, params, testState, st)
	utils.ResetOps()
	var out StateHE
	for b := 0; b < 16; b++ {
		if out[b], err = ae.SubByteV2(st[b]); err != nil {
			t.Fatalf("SubByteV2 %d: %v", b, err)
		}
	}
	ops := utils.Ops
	t.Logf("SubBytesV2 (16 bytes): %d relin, %d rescale (%d relin, %d rescale per byte)",
		ops.Relin, ops.Rescale, ops.Relin/16, ops.Rescale/16)
	// Non-lazy BCKK baseline: exactly one relin + one rescale per degree-2+ monomial,
	// i.e. 255 - 8 = 247 per byte (no alignment squaring).
	if ops.Relin != 16*247 || ops.Rescale != 16*247 {
		t.Errorf("op count = %d relin / %d rescale, want %d / %d (247 per byte)",
			ops.Relin, ops.Rescale, 16*247, 16*247)
	}
	dbgState("out", ae, params, want, out)

	checkState(t, params, ecd, dec, out, want, "SubBytesV2")
}

func TestSubBytesV3(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)

	st, err := EncStateRepl(params, ecd, enc, testState, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc state: %v", err)
	}

	want := testState
	SubBytes(want[:])

	dbgState("in", ae, params, testState, st)
	utils.ResetOps()
	var out StateHE
	for b := 0; b < 16; b++ {
		if out[b], err = ae.SubByteV3(st[b]); err != nil {
			t.Fatalf("SubByteV3 %d: %v", b, err)
		}
	}
	ops := utils.Ops
	t.Logf("SubBytesV3 (16 bytes, lazy): %d relin, %d rescale (%d relin, %d rescale per byte)",
		ops.Relin, ops.Rescale, ops.Relin/16, ops.Rescale/16)
	// Lazy baseline: 61 reused-factor monomials + 8 output-bit accumulators = 69 per byte.
	if ops.Relin != 16*69 || ops.Rescale != 16*69 {
		t.Errorf("op count = %d relin / %d rescale, want %d / %d (69 per byte)",
			ops.Relin, ops.Rescale, 16*69, 16*69)
	}
	dbgState("out", ae, params, want, out)

	checkState(t, params, ecd, dec, out, want, "SubBytesV3")
}

func TestSubBytesV4(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)

	st, err := EncStateRepl(params, ecd, enc, testState, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc state: %v", err)
	}

	want := testState
	SubBytes(want[:])

	dbgState("in", ae, params, testState, st)
	utils.ResetOps()
	var out StateHE
	for b := 0; b < 16; b++ {
		if out[b], err = ae.SubByteV4(st[b]); err != nil {
			t.Fatalf("SubByteV4 %d: %v", b, err)
		}
	}
	ops := utils.Ops
	t.Logf("SubBytesV4 (16 bytes, lazy deg>=4): %d relin, %d rescale (%d relin, %d rescale per byte)",
		ops.Relin, ops.Rescale, ops.Relin/16, ops.Rescale/16)
	// 61 factors + 29 low-degree (2,3) leaves = 90 build, + 8 output-bit accumulators = 98 per byte.
	if ops.Relin != 16*98 || ops.Rescale != 16*98 {
		t.Errorf("op count = %d relin / %d rescale, want %d / %d (98 per byte)",
			ops.Relin, ops.Rescale, 16*98, 16*98)
	}
	dbgState("out", ae, params, want, out)

	checkState(t, params, ecd, dec, out, want, "SubBytesV4")
}

func TestShiftRows(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)

	st, err := EncStateRepl(params, ecd, enc, testState, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc state: %v", err)
	}

	want := testState
	ShiftRows(want[:])

	dbgState("in", ae, params, testState, st)
	out := ShiftRowsHE(st)
	dbgState("out", ae, params, want, out)

	checkState(t, params, ecd, dec, out, want, "ShiftRows")
}

func TestMixColumnsV2(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)

	st, err := EncStateRepl(params, ecd, enc, testState, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc state: %v", err)
	}

	want := testState
	MixColumns(want[:])

	dbgState("in", ae, params, testState, st)
	out, err := ae.MixColumnsV2(st)
	if err != nil {
		t.Fatalf("MixColumnsV2: %v", err)
	}
	dbgState("out", ae, params, want, out)

	checkState(t, params, ecd, dec, out, want, "MixColumnsV2")
}

func TestMixColumnsV1(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)

	st, err := EncStateRepl(params, ecd, enc, testState, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc state: %v", err)
	}

	want := testState
	MixColumns(want[:])

	dbgState("in", ae, params, testState, st)
	out, err := ae.mixColumnsV1(st)
	if err != nil {
		t.Fatalf("mixColumnsV1: %v", err)
	}
	dbgState("out", ae, params, want, out)

	checkState(t, params, ecd, dec, out, want, "MixColumnsV1")
}
