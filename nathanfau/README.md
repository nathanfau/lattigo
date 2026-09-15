# TODO 

LastRoundV2
EvalMod in 7lvl (using the cosine's parity)
ModSwitch
FreeXor (for less than 2^15 aes blocks, new packing, new BTS...)
rounds without oracle call
logN 16

# Lattigo Fork

This repository is a fork of Lattigo v6 made on the 4th of July 2026.

The goal of the code is homomorphic **AES-128 transciphering** under CKKS: an encrypted AES state is advanced round by round, and the level budget consumed by each round is given back by an `IntRootBoot` refresh (Algorithm 1 of [BCKK25]).

## Code overview

Listed base first: each package only depends on the ones above it.

- **`utils`** — Generic HE helpers: level alignment, balanced reduction trees, bit-distance
  precision measurement, box-drawing result tables, timing stats.
- **`params2`** — The CKKS and bootstrapping parameter set of the transciphering pipeline.
  Every other package builds its own parameters in its test.
- **`debug`** — Printing and inspection: slot dumps (CI and Std), modulus-chain traces, precision tables over one or many ciphertexts, parameter summaries.
- **`cleaning`** — Pushes noisy near-bit slots back towards 0/1: three strengths.
- **`trigo`** — Homomorphic cos, sin and exp by Chebyshev approximation.
- **`bitbatching`** — Packs several bit planes into one ciphertext and extracts them back out.
- **`convctx`** — Conjugate-invariant ↔ standard domain switch at full capacity ([BCKK25], Fig. 1 and 2).
- **`blockpack`** — Packs AES blocks into 8×8 ciphertexts, 8 groups of 8 bits; each group folds two bytes of the state into the two halves of the slot vector.
- **`aes`** — Bit-sliced AES-128 under CKKS: SubBytes in three variants, ShiftRows, MixColumns,
  AddRoundKey, KeyExpansion. Ships a cleartext AES next to it, used as the test oracle.
- **`algo1`** — Algorithm 1 of [BCKK25], the conjugate-invariant `IntRootBoot`.
- **`bbbts`** — Algorithm 1 of [BKSS24], the "slim"-order bootstrap whose `EvalMod` is replaced by a direct evaluation of the root of unity.
- **`transciphering`** — Wires the `algo1` refresh into the AES round pipeline.

### The pipeline

`transciphering` is where everything meets. One AES round is

```
SubBytes → Refresh (algo1, ShiftRows applied during the pause) → MixColumns → AddRoundKey → Cleaning
```

A full AES-128 is `FirstRound` + 9 × `Round` + `LastRoundV1`, which is what `TestAES` runs.

## Documentation

*(to be written)*

## Differences with the doc

*(to be written)*

### Discussion: XOR and cleaning, average vs worst case (see tests in the cleaning package)

Squaring a noisy bit cleans it exactly half the time. Take `(x-y)^2` on `x = bx + ex`,
`y = by + ey`. When the two bits agree the integer parts cancel and the result is
`(ex-ey)^2`, an error squared: free cleaning. When they differ the result is
`1 +/- 2(ex-ey)`, so the error becomes `2|ex-ey|`: twice the input error typically,
four times it when `ex` and `ey` reinforce instead of cancelling.

Half the slots therefore come out cleaned and half come out amplified. Over thousands
of slots the extreme of that second half is always drawn, so the worst slot pays the
full `4e`. The worst slot is the only thing that decides whether
a transciphering round fails.

That is the whole tension: the square buys a lot of average precision from the half
that cleans itself, and pays for it exactly where it hurts. The average even lands
where the split predicts: at input 2^-12 one half of the slots goes to ~24 bits and
the other stays near 11, for a mean around 18.

| raw XOR, input 2^-12 | avg | worst |
|---|---|---|
| `XorSq` `(x-y)^2` | 18.56 | 9.56 |
| `XorNoSq` `x+y-2xy` | 12.78 | 10.55 |

`XorNoSq` cleans nothing: whatever the two bits, the result carries `|ex +/- ey|`, at
most `2e`. No free average, but -1 bit per gate instead of -2.

**With XorNoSq.** For `y` an exact bit (ie fresh ct), `h(x XOR y) = h(x) XOR y`: the cleaning
polynomials are smoothstep, so `h(1-t) = 1-h(t)`, and `XorNoSq(x,y)` is `x` or
`1-x`. Cleaning before, cleaning after and cleaning only `x` therefore give the same
digits. Since the second operand of `AddRoundKey` is a fresh round key, the existing
`AddRoundKey` then `Clean` is already optimal, and it is also the cheapest
(~44 ms against ~85 ms for cleaning both operands). Writing the polynomial out by
hand instead of calling `Cleaning` is worth +0.05..0.14 bit of average and 2-7% of
time (we do not keep this code path for the moment).

**If XorSq is kept instead.** The identity breaks, and every choice starts to
matter. Cleaning before beats cleaning after by the square's own worst-slot cost:
+2 bits on two noisy operands, +1 bit against a round key. Cleaning only `x` then
matches cleaning both for half the polynomial work. And the basin becomes a real
hazard: at input 2^-2, `XorSq` followed by cleaning reaches worst = -3.35 bits
(absolute error ~10) because the XOR pushed values out of [0,1] and the polynomial
diverged instead of contracting; cleaning first keeps both operands inside it.

### `convctx`: masking and key switch are swapped, on purpose

[BCKK25] Fig. 1 masks first (step 2) and key switches after (step 3). We do the opposite:
lattigo's `RealToComplex` fuses the ring embedding of step (1) with the key switch of step (3)
into a single call, so taking the paper's order means unfolding by hand and applying the
ring-swap key separately. The two orders give the same result but not the same noise, because 
the key switch then runs one level higher, with one gadget decomposition digit more.

Both orders were implemented and measured against the exact value, 16 draws, at the two levels
where the pipeline actually converts CI → Std (2, the nibbles out of SubBytes; 21, out of
`EvalCos`/`EvalSin` in algo1). Only ours survives in the tree, as `CIToStandard`; the paper-order
variant was deleted once the table below settled the question:

| chain | level | ours, avg / worst | paper, avg / worst | paper − ours |
|---|---|---|---|---|
| params2, logN 11 | 2  | 29.25 / 26.60 | 28.94 / 26.47 | −0.31 / −0.13 |
| params2, logN 11 | 21 | 29.25 / 26.46 | 28.94 / 26.53 | −0.31 / +0.07 |
| params2, logN 12 | 2  | 28.74 / 26.10 | 28.44 / 25.86 | −0.30 / −0.24 |
| params2, logN 12 | 21 | 28.74 / 26.00 | 28.44 / 25.90 | −0.30 / −0.10 |

Our order wins ~0.30 bit on the average everywhere, and on the worst slot in three cases out of four; the remaining `+0.07` sits inside the scatter
of a worst-slot statistic. So the swap is kept.

### `dnum` is not a free parameter here, so `P` is sized in `Q` primes

The number of key-switching digits never appears as a parameter in lattigo.
It is derived — `dnum = Ceil(#Q/#P)` (`core/rlwe/params.go:543`) — and the decomposition itself is
built on `nbPi := len(P)` **consecutive `Q` primes per digit** (`ring/basis_extension.go:336`), the
last digit taking the remainder. The digit width is therefore counted in *primes*, and the only way
to set `dnum` is to choose how many primes `P` holds.

So `dnum` cannot be set by hand here. Parameter tables in [BCKK25] are written the other way
round: `dnum` is the input, the chain is cut into `dnum` digits, and `P` is then sized *in bits* to
cover the widest one. OpenFHE exposes that directly (`numLargeDigits`), and their Table 1
shows HEaaN does too — its `Param` rows use **62-bit** `P` primes against a 33/38-bit chain and
print `dnum = 3` with 6 `P` primes over ~32 `Q` primes, where `Ceil(32/6) = 6`. The two conventions
coincide **iff the `P` primes have the size of the `Q` primes** — only then is "as many primes as
`P` has" the same cut as "as many bits as `P` has".

Their XBOOT rows are exactly that case, and there lattigo's formula is exact: 38-bit `P` primes on a
36–38-bit chain give (17 `Q` primes, 5 `P`) → 4, (18, 4) → 5, (20, 2) → 10, the three `dnum`
printed, with `log(PQ)` adding up to the bit (829, 828).

**Hence our `P` departs from theirs**, and it has to: lattigo caps a `P` prime at 61 bits
(`checkModuliLogSize`, `MaxModuliSize+1`), so `62*6` is refused outright, and even at 61 bits 6
primes would mean `dnum = 6`, not 3. The faithful translation preserves `log(P)` and the digit,
never the prime count.

We follow XBOOT instead and size `P` in `Q` primes, which lands `log(P)` just above the digit by
construction — which is where it belongs, because the margin is a step and not a slope: swept at
`logN = 12`, every `P` above its digit returns the same precision to the tenth of a bit, and below
it the worst slot falls by about a bit per bit. Every bit past the digit is `log(QP)` given away.
Our old `6*61+50` sat at +124; the current `8*42` sits at +10.

The rule that follows: pick `dnum`, then take `#P = Ceil(#Q/dnum)` primes of `Ceil(digit/#P)` bits.
More `P` primes at that `dnum` pay `#Q + #P` for nothing; fewer raise `dnum`. On our chain,
`dnum 5 → 7*42`, `dnum 4 → 8*42`, `dnum 3 → 11*41`.

What we keep from [BCKK25] is the chain `Q` — every level, every scale, every calibration of the
pipeline is theirs. `P` carries no message and no level, so its shape is invisible to the circuit
and changes only what a key-switch costs and how much noise it leaves. There is no faithfulness to
lose here, only bits.

## Testing

Every package is tested individually and builds its own parameters, so tests can be run in any
order:

```bash
go test ./nathanfau/<package>/ -v
```

`debug`, `params2` and `utils` have no test of their own; they are exercised through the packages that use them. Three packages take flags.

### algo1

`TestAlgo1Zone` is `TestAlgo1` on the transciphering chain with a zone: the primes from
Conv_{Real->Cplx} up to the bit extraction, and the scale they run at, take `-zone` bits (default
33; 38 is the chain without a zone). The scale is checked at each edge of the zone. `-pre` and
`-post` add that many XOR layers before the refresh and after the bit extraction (default 0), so
the scale drifts the way the AES circuit makes it drift.

```bash
go test ./nathanfau/algo1/ -run '^TestAlgo1Zone$' -v -zone 33 -pre 3 -post 4 -timeout 0
```

### blockpack

`-logn` sets the CI ring log-degree, hence the capacity (default 12, so 2^11 blocks).
This layout is built for 2^15 AES blocks, that is `-logn 16`. A variant for fewer blocks will be
available on this repo soon.

```bash
go test ./nathanfau/blockpack/ -v
```

### transciphering

`-rounds` is how many AES middle rounds `TestTransciphering` chains (default 1), `-subbytes` picks the SubBytes variant, 1 to 3 by decreasing cost (247, 98 and 69 relinearisations per byte). The default is 2, the middle one.

Five more flags select the variants the pipeline runs on, so a whole matrix of circuits can be
compared without touching the code:

| flag | values | default |
|---|---|---|
| `-xor` | `nosq` for `x+y-2xy`, `sq` for `(x-y)^2` | `nosq` |
| `-clean` | `cleaning` (2 levels), `smoother` or `verysmoother` (3 levels, one prime more) | `cleaning` |
| `-place` | `after` (AddRoundKey then Cleaning), `both` (clean both operands, then XOR), `one` (clean the state only, then XOR) | `after` |
| `-seed` | seed of the block draw, 0 draws one from the clock | `0` |
| `-zone` | size of the primes of the AES circuit and the bit extraction, and of their scale; `0` or `38` is no zone | `0` |

`-xor` is the one that decides whether the cipher completes. It applies to **every** gate: the
initial AddRoundKey, the round AddRoundKey, and the XOR trees of MixColumns. `sq` fails a full
AES-128 at every combination of the other flags, because the extra bit it costs on the worst slot
accumulates until the cleaning polynomial leaves its basin and amplifies instead of contracting;
`nosq` holds a stationary error over the ten rounds, which is why it is the default.

`-clean` is baked into the parameters rather than switched later: the polynomial's depth decides
how many primes the chain carries. `-place` does not change the level map, but it does change the
level a round key has to be encrypted at, which `Context.ARKKeyLv` and `Context.LastKeyLv` carry.

`-seed` fixes which blocks the batch carries, not the key or the encryption noise, which lattigo
draws from `crypto/rand`: two runs on one seed are far closer than on two, not identical.

`-zone`, like `-clean`, is baked into the parameters. With a zone, the state and the round keys
are encrypted at its scale, the refresh leaves it for the bootstrap and Algo1 lands back on it;
without one the pipeline runs exactly as before. A run in a zone shows in the CSV through
`q_sizes` and `chain_id`.

```bash
go test ./nathanfau/transciphering/ -run '^TestTransciphering$' -v -subbytes 3 -rounds 1 -timeout 0
go test ./nathanfau/transciphering/ -run '^TestAES$' -v -subbytes 2 -timeout 0
go test ./nathanfau/transciphering/ -run '^TestAES$' -v -xor sq -clean verysmoother -place both -seed 1 -timeout 0
```

`-timeout 0` matters as soon as the run gets long: a whole transciphering exceeds Go's default 10-minute budget.
`TestAES` runs the 10 rounds of a full AES-128, takes several minutes, and is skipped under `-short`.
Both check every block of the batch against the cleartext AES oracle after every step, so `-v` prints a per-step precision and correctness trace.

Note that key generation, not the circuit, is what dominates memory: on a small machine, lower
`logN` before lowering anything else.

---

LLMs were used occasionally to help write or refactor code, and to assist with the documentation