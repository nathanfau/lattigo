package aes

// We deliberately do not test:
//-xor: it has no test of its own because MixColumns already covers it.
//-ShiftRows: ShiftRows consumes the nibbles ALREADY split into their
//  Re / Im streams by the Algo1 pause, so exercising it here would mean
//  reproducing that split outside the refresh it belongs to.

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/nathanfau/blockpack"
	"github.com/tuneinsight/lattigo/v6/nathanfau/debug"
	"github.com/tuneinsight/lattigo/v6/nathanfau/utils"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// testBlocks fills the whole batch with independent random AES blocks, one per slot,
//
//	and logs the seed so a failure can be replayed.
func testBlocks(t *testing.T, params ckks.Parameters) [][16]byte {
	t.Helper()
	seed := time.Now().UnixNano()
	t.Logf("random seed = %d", seed)
	return blockpack.RandomBlocks(params, rand.New(rand.NewSource(seed)))
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
	debug.DbgParams("aes", params)
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

// worstPrec is the precision of the worst slot of the whole packed state, in bits.
func worstPrec(t *testing.T, ecd *ckks.Encoder, dec *rlwe.Decryptor, p blockpack.Packed, want [][]float64) float64 {
	t.Helper()
	st, err := utils.BitDistance(ecd, dec, p.Cts(), want, 0)
	if err != nil {
		t.Fatalf("BitDistance: %v", err)
	}
	agg, _ := utils.WorstOf(st)
	return -math.Log2(agg.Worst)
}

// subByteVariant is one SubBytes implementation under test, with its expected
// relin==rescale count per S-box circuit (the BCKK op-count baseline).
type subByteVariant struct {
	name        string
	fn          func(ByteHE) (ByteHE, error)
	wantPerByte int
}

// TestSubBytes runs every SubBytes variant on the same batch of blocks as sub-tests
func TestSubBytes(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)

	blocks := testBlocks(t, params)
	st, err := blockpack.Encrypt(params, ecd, enc, blocks, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc state: %v", err)
	}
	want := make([][16]byte, len(blocks))
	for i, blk := range blocks {
		want[i] = blk
		SubBytes(want[i][:])
	}

	variants := []subByteVariant{
		{"V2", ae.SubByteV2, 247}, // non-lazy: one relin+rescale per degree-2+ monomial (255-8)
		{"V3", ae.SubByteV3, 69},  // lazy: 61 reused factors + 8 output-bit accumulators
		{"V4", ae.SubByteV4, 98},  // lazy deg>=4: 61 factors + 29 low-degree leaves + 8 accumulators
	}

	var leads []string
	var rows [][]string
	for _, v := range variants {
		t.Run(v.name, func(t *testing.T) {
			utils.ResetOps()
			var out blockpack.Packed
			t0 := time.Now()
			for g := range 8 {
				var e error
				if out[g], e = v.fn(st[g]); e != nil {
					t.Fatalf("%s group %d: %v", v.name, g, e)
				}
			}
			dur := time.Since(t0)
			ops := utils.Ops

			if e := blockpack.Check(params, ecd, dec, out, want); e != nil {
				t.Fatalf("%s %v", v.name, e)
			}
			if ops.Relin != 8*v.wantPerByte || ops.Rescale != 8*v.wantPerByte {
				t.Errorf("%s op count = %d relin / %d rescale, want %d / %d (%d per circuit)",
					v.name, ops.Relin, ops.Rescale, 8*v.wantPerByte, 8*v.wantPerByte, v.wantPerByte)
			}

			leads = append(leads, v.name)
			rows = append(rows, []string{
				fmt.Sprintf("%d", ops.Relin/8), fmt.Sprintf("%d", ops.Relin),
				fmt.Sprintf("%d", ops.Rescale/8), fmt.Sprintf("%d", ops.Rescale),
				dur.Round(time.Millisecond).String(),
				fmt.Sprintf("%.2f", worstPrec(t, ecd, dec, out, blockpack.WantSlots(params, want))),
			})
		})
	}

	groups := []utils.TableGroup{
		{Name: "relin", Cols: []string{"/circ", "total"}},
		{Name: "rescale", Cols: []string{"/circ", "total"}},
		{Name: "circuit", Cols: []string{"time", "worst prec(b)"}},
	}
	fmt.Printf("SubBytes variants\n%s", utils.BoxTable("SubBytes", groups, leads, rows))
}

func TestMixColumns(t *testing.T) {
	params, ecd, enc, dec, ae := newAESContext(t)

	blocks := testBlocks(t, params)
	st, err := blockpack.Encrypt(params, ecd, enc, blocks, params.MaxLevel())
	if err != nil {
		t.Fatalf("enc state: %v", err)
	}

	want := make([][16]byte, len(blocks))
	for i, blk := range blocks {
		want[i] = MixColumnsRM(blk)
	}

	debug.DbgBitStd("in bit[0][0]", ae.eval, st[0][0], blockpack.SlotVec(params, blocks, 0, 0))

	out, err := ae.MixColumns(st)
	if err != nil {
		t.Fatalf("MixColumns: %v", err)
	}

	debug.DbgBitStd("out bit[0][0]", ae.eval, out[0][0], blockpack.SlotVec(params, want, 0, 0))

	if err := blockpack.Check(params, ecd, dec, out, want); err != nil {
		t.Fatalf("MixColumns %v", err)
	}
}
