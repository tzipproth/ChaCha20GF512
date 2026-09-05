# ChaCha20GF512
## Maximizing What Can Be Proven About a PRNG

> **ChaCha20 is already an excellent computational PRNG. What exact mathematical properties can be added on top of it?**

ChaCha20GF512 starts with ChaCha20 and adds a second, independent notion of randomness: **exact 512-wise independence of 64-bit output words**.

The goal is not to fix a statistical weakness in ChaCha20. No such weakness is assumed. The point is that computational pseudorandomness and exact finite-order independence are different properties, and they can coexist in one generator.

For a fixed public stream identifier, ordinary 256-bit-key ChaCha20 has at most `2^256` keyed streams. That is already enough for cryptographic pseudorandomness, but it is too small a family to make even five selected 64-bit outputs *exactly* jointly uniform over all `2^320` possible 5-tuples.

ChaCha20GF512 adds an independently seeded degree-511 polynomial over `GF(2^128)`. Each 128-bit field evaluation supplies two adjacent 64-bit GF words (low half, then high half), which are XORed with the corresponding ChaCha20 words. The resulting family has the exact guarantee:

> **For any 512 distinct `uint64_t` output positions, the 512 words are jointly uniform when the GF coefficients are sampled uniformly and independently.**

This is not a statistical-test statement. It follows from finite-field interpolation.

### What ChaCha20GF512 adds to ChaCha20

| Property | ChaCha20 | ChaCha20GF512 |
|---|---|---|
| Computational pseudorandomness | Yes, this is ChaCha20's main role | Retained as the cryptographic component |
| Exact 512-wise independence of 64-bit words | No such exact guarantee follows from the 256-bit keyed family | **Yes** |
| Exact theorem valid at arbitrary selected positions | Not an information-theoretic claim of ChaCha20 | **Yes, for any <= 512 distinct positions** |
| Addressable words in one logical stream | Original ChaCha: `2^64` blocks = `2^67` `uint64_t` words per fixed nonce | **`2^128` `uint64_t` words** |
| Random access | Yes, counter based | Yes |
| Parallel evaluation of one logical stream | Yes in principle | Yes, implemented on 1024-output-word / 512-evaluation FFT boundaries |

RFC 8439 notes that the original ChaCha design used a 64-bit nonce and a 64-bit block counter. With eight 64-bit words per ChaCha block, that gives `2^67` 64-bit words per fixed nonce. ChaCha20GF512 instead uses HChaCha20 for stream separation and exposes a 128-bit **word-position** space, supporting exactly `2^128` 64-bit word positions in one logical stream. The ChaCha block index is the word position divided by eight and is stored in state words 12..15.

The supported stream is therefore longer by a factor of

```text
2^128 / 2^67 = 2^61
```

than the original 64-bit-counter ChaCha layout for one fixed nonce.

The implementation does not silently wrap at the end of that position space.

---

## Is this “more random” than ChaCha20?

There is no single scalar definition that totally orders good PRNGs as “more random” or “less random.”

Different notions measure different things:

- **computational pseudorandomness:** can an efficient algorithm distinguish the stream from random?
- **entropy / seed-space size:** how many initialized streams are possible?
- **exact k-wise independence:** are arbitrary sets of up to `k` outputs *exactly* jointly uniform over the initialization family?
- **period / addressable stream length:** how long can the generator run before its position space is exhausted?
- **statistical-test behavior:** does a finite sample trigger a particular empirical test?

ChaCha20GF512 does **not** claim a known improvement over ChaCha20 in computational pseudorandomness. ChaCha20 is already used precisely because its output is computationally pseudorandom.

But in exact finite-order independence the hybrid is genuinely stronger:

> **ChaCha20GF512 adds an exact information-theoretic distribution property that a 256-bit-key ChaCha20 family with fixed/public stream ID cannot possess.**

So “more random” is too vague, but **“strictly stronger in this exact mathematical sense”** is accurate.

The construction is

```text
j = floor(i / 2)
G(2j)   = low64(P(j))
G(2j+1) = high64(P(j))
R(i)    = ChaCha20(K, stream_id, i) XOR G(i)
```

where

```text
P(x) = a0 + a1*x + ... + a511*x^511
```

is a random polynomial over `GF(2^128)`.

The 512 GF coefficients occupy 8192 bytes. Together with the 32-byte ChaCha key, the full explicit seed is 8224 bytes.

---

## Stream length / period

ChaCha20GF512 is best viewed as a finite, position-indexed stream rather than as a small-state recurrence with a conventional proven minimal period.

Valid word positions are

```text
0 <= i < 2^128.
```

That is

```text
2^128 uint64_t words
= 2^131 bytes
~= 2.72e39 bytes.
```

Even at an intentionally absurd sustained output rate of **1 exabyte per second**, exhausting the position space would take about

```text
8.6e13 years.
```

The implementation treats the end of this space as exhaustion:

- `next_int()` does not silently wrap to position zero;
- random access uses a 128-bit position;
- parallel bulk generation validates and splits ranges within the same 128-bit domain.

This should not be confused with a theorem that the **minimal period of the 64-bit output values** is exactly `2^128`. No such claim is needed. Equal 64-bit output values may occur at different positions, just as they do in genuine random data.

The useful guarantee is simpler:

> **There are `2^128` distinct supported positions, no position is reused, and the exact 512-wise theorem applies over that entire domain.**

---

# 1. Construction

## 1.1 128-bit position

The implementation uses the small POD type

```cpp
ChaCha20GF512Position128
```

with two 64-bit limbs:

```text
Position = hi * 2^64 + lo.
```

The old 64-bit APIs remain available as convenience overloads, so existing code using positions below `2^64` does not need to change.

The 128-bit output-word position serves three roles:

1. divided by two, it selects the `GF(2^128)` polynomial evaluation point;
2. its low bit selects the low or high 64-bit half of that field value;
3. divided by eight, it gives the ChaCha block index because each ChaCha block contains eight `uint64_t` words.

---

## 1.2 ChaCha20 component

The implementation uses the 20-round ChaCha core with a 256-bit key.

RFC 8439 standardizes the later IETF layout with a 32-bit block counter and 96-bit nonce, while explicitly noting that the original ChaCha design used a **64-bit nonce and 64-bit block counter**.

The current ChaCha20GF512 revision goes one step further for its position-indexed use case:

1. the public 64-bit `stream_id` is placed in a 128-bit HChaCha20 input (`stream_id || 0`);
2. HChaCha20 derives a 256-bit subkey;
3. all four ChaCha state words 12..15 are then used as one 128-bit block index.

HChaCha20 uses the standard 20 ChaCha rounds without feed-forward and extracts words

```text
0, 1, 2, 3, 12, 13, 14, 15.
```

This keeps independent logical streams through subkey derivation while leaving all four counter words available. The current `2^128`-word output domain uses `2^125` ChaCha blocks, so it remains well inside that counter capacity.

This is a **project-specific position-indexed construction**, not the RFC 8439 IETF layout and not a claim that the whole generator is a standardized XChaCha mode. HChaCha20 itself is the established subkey primitive used by XChaCha constructions.

ChaCha20GF512 uses ChaCha as a counter-indexed pseudorandom component, not as a state-updating DRBG. There is no periodic rekeying. Random access and deterministic parallel evaluation are deliberate design goals.

---

## 1.3 Galois field GF(2^128) component

The second component is a degree-at-most-511 polynomial

```text
P(x) = a0 + a1*x + ... + a511*x^511
```

over `GF(2^128)`.

There are 512 independent coefficients:

```text
512 * 128 bits = 65536 bits = 8192 bytes.
```

The implementation uses the binary field defined by

```text
x^128 + x^7 + x^2 + x + 1,
```

the polynomial used for the `GF(2^128)` multiplication in GHASH/GCM.

With the implementation's low-bit polynomial representation, reduction of an overflowed `x^128` term uses the constant

```text
0x87
```

because

```text
x^128 = x^7 + x^2 + x + 1.
```

Field addition is bitwise XOR. Field multiplication is carry-less polynomial multiplication followed by reduction modulo the degree-128 polynomial.

On x86/x86-64 the implementation uses `PCLMULQDQ` when available; otherwise it uses a portable shift/XOR implementation.

---

## 1.4 Two 64-bit output words per GF evaluation

For output-word position `i`, define

```text
j = floor(i / 2).
```

The mathematical field value is

```text
P(j) in GF(2^128).
```

Both halves are used:

```text
G(2j)   = low64(P(j))
G(2j+1) = high64(P(j)).
```

So no half of a field evaluation is discarded. A uniform `GF(2^128)` value is exactly a uniform 128-bit vector, and its low and high halves are jointly uniform 64-bit words. Selecting only one half is a surjective linear projection; selecting both halves is the identity on the 128 output bits.

This mapping is also the main performance improvement over the previous revision: one 512-point additive FFT now supplies 1024 `uint64_t` GF words.

---

## 1.5 Combination

The final output is

```text
R(i) = C(i) XOR G(i),
```

where

- `C(i)` is the 64-bit ChaCha word at output position `i`;
- `G(2j) = low64(P(j))`;
- `G(2j+1) = high64(P(j))`.

For any fixed ChaCha stream, XOR is simply a position-dependent translation of the GF-derived words. A translation is a bijection, so it cannot destroy the exact finite-order uniformity supplied by the GF component.

The intended division of labor is:

| Component | Role |
|---|---|
| ChaCha20/HChaCha20 | computational pseudorandomness |
| degree-511 GF polynomial | exact finite-order distribution theorem |
| low/high half mapping | turns each uniform 128-bit evaluation into two adjacent 64-bit words without discarding either half |
| XOR | lets both properties coexist in one output stream |

---

# 2. What is proven exactly

## 2.1 Full GF evaluations are exactly 512-wise independent

Let

```text
P(x) = a0 + a1*x + ... + a511*x^511
```

where all 512 coefficients are independent and uniformly distributed in `GF(2^128)`.

For any `t <= 512` distinct field points

```text
x1, x2, ..., xt,
```

the vector

```text
(P(x1), P(x2), ..., P(xt))
```

is exactly uniform over

```text
GF(2^128)^t.
```

This is the standard random-polynomial construction of `k`-wise independent variables.

For `t = 512`, prescribing arbitrary field values at 512 distinct points determines exactly one degree-at-most-511 polynomial by interpolation.

For smaller `t`, the remaining coefficient freedom gives the same uniform result.

This is an exact algebraic statement, not a statistical-test result.

---

## 2.2 The 64-bit output words inherit exact 512-wise independence

Take any `t <= 512` distinct output-word positions

```text
i1, i2, ..., it.
```

Map them to GF evaluation points

```text
jr = floor(ir / 2).
```

Several selected words may refer to the same GF point, but at most two can do so: one asks for the low half and one for the high half. Let there be `r <= t <= 512` distinct GF points among the selected positions.

By Section 2.1, the `r` full values

```text
(P(j1), ..., P(jr))
```

are independent and uniform 128-bit field elements. For each such value, the requested output coordinates are either its low 64 bits, its high 64 bits, or both. Coordinate selection from a uniform 128-bit vector is uniform on the selected 64 or 128 bits.

Therefore

```text
(G(i1), ..., G(it))
```

is exactly uniform over

```text
({0,1}^64)^t.
```

For a fixed ChaCha stream,

```text
(R(i1), ..., R(it))
```

is exactly uniform as well, because XOR with the fixed ChaCha vector is a bijection.

So the production theorem remains:

> **Any selection of at most 512 distinct `uint64_t` positions in the supported `2^128`-word stream is jointly uniform over the full GF coefficient seed space.**

---

## 2.3 Exact seed-space counting

The GF seed contains

```text
65536 bits.
```

For `t <= 512` selected 64-bit output positions, every possible `t`-tuple has exactly

```text
2^(65536 - 64*t)
```

GF-seed preimages, for every fixed ChaCha stream.

In particular, for exactly 512 selected words, every possible 32768-bit tuple occurs exactly `2^32768` times over the complete GF seed space.

There is also a useful stronger **structured** statement. Take 512 distinct GF evaluation points and expose both halves of each. The resulting 1024-word / 65536-bit tuple is exactly uniform, and every such tuple has exactly one GF-seed preimage. For example, output positions `0..1023` have this property.

This does **not** claim arbitrary 1024-wise independence; it says that the full 128 bits of each of 512 independent field evaluations are no longer discarded.

---

## 2.4 Arbitrary positions, not only consecutive outputs

The positions do not need to be adjacent.

For example, any set such as

```text
0
17
4711
2^64 + 123
2^100 + 7
...
```

is valid as long as all selected output positions are distinct and below `2^128`.

Unlike the previous low-half-only revision, two adjacent word positions now intentionally share one interpolation point: positions `2j` and `2j+1` use the low and high halves of `P(j)`. This does not weaken the theorem because a uniform 128-bit evaluation provides two jointly uniform 64-bit coordinates.

Across the complete supported `2^128`-word stream, the GF component uses evaluation points `0 <= j < 2^127`.

---

## 2.5 Seed-space cost and information use

Exact 512-wise independence of arbitrary 64-bit words requires at least

```text
512 * 64 = 32768 seed bits
```

in any exact family.

ChaCha20GF512 uses

```text
512 * 128 = 65536 GF seed bits.
```

So the seed is still larger than the information-theoretic minimum needed for the **arbitrary-512-word** theorem. The 128-bit coefficients are retained because the construction works over `GF(2^128)`, which gives a huge evaluation domain and maps naturally to PCLMUL-accelerated arithmetic.

The previous revision then threw away half of every 128-bit evaluation. The current revision does not: 512 evaluations provide 1024 output words, so an aligned 1024-word block can expose all 65536 bits of GF-seed dimension exactly.

At 8192 bytes, the GF coefficient vector remains small in ordinary software terms.

---

# 3. Where the exact guarantee stops

The full `GF(2^128)` polynomial evaluations are exactly 512-wise independent as 128-bit variables and are not a general cryptographic PRG. Once enough full evaluations are known, the polynomial is algebraically reconstructible.

For 513 selected **full 128-bit field evaluations**, there is necessarily a non-trivial GF-linear relation because the polynomial has only 512 field coefficients.

The production stream now exposes both halves of an evaluation at adjacent positions. Consequently, 1024 appropriately paired words can reveal 512 complete evaluations, and an additional complete evaluation makes the GF-only algebraic dependence explicit. For example, 1026 consecutive words starting at an even position contain 513 complete field evaluations.

That still does **not** prove that arbitrary 513 distinct 64-bit output positions fail to be independent. A set of 513 words can choose only one 64-bit coordinate from each of 513 distinct field points, and projection can hide full-field linear relations. Therefore this README keeps the conservative theorem that is directly proved:

> **At least 512-wise independence of arbitrary 64-bit output words.**

The GF seed contains 65536 bits, so pure seed counting rules out exact arbitrary independence beyond 1024 64-bit words. The exact maximum order between 512 and 1024 remains a separate algebraic question.

Beyond the finite exact guarantee, the intended pseudorandomness argument comes from the ChaCha/HChaCha component.

---

# 4. What relies on the ChaCha/HChaCha assumption

The exact 512-wise result is unconditional once the GF coefficients are sampled as specified.

The computational claim is different.

The current construction assumes the 20-round HChaCha20 subkey derivation and the 20-round ChaCha core behave as secure pseudorandom primitives in this position-indexed use.

This is not a proof of ChaCha20 or HChaCha20 security, and this exact 128-bit-counter composition is not itself a standardized primitive.

What XOR gives us is a clean robustness argument once the ChaCha component is modeled as pseudorandom and the GF seed is independent:

- if the ChaCha stream is computationally indistinguishable from random, XOR with an independently generated GF stream cannot make it efficiently distinguishable merely by introducing a fixed reversible translation in the reduction;
- if one conditions on any fixed ChaCha stream, the exact GF 512-wise theorem remains true.

The two guarantees therefore have different failure modes.

### Graceful degradation

- A future cryptanalytic break of the ChaCha/HChaCha assumption would not invalidate the exact 512-wise GF theorem.
- Ignoring the GF theorem beyond its proven order does not remove the computational pseudorandom component.

This does **not** mean the hybrid would necessarily remain cryptographically secure after a break of ChaCha/HChaCha. It means only that the information-theoretic bounded-independence statement is structurally separate.

---

# 5. The guarantee is over a generator family

A fixed-seed PRNG instance is deterministic.

Therefore

> “ChaCha20GF512 is 512-wise independent”

is shorthand for a statement about the family obtained by sampling the GF coefficient vector uniformly.

Equivalently, it is an exact counting statement over the complete GF seed space.

For every fixed ChaCha key and fixed public `stream_id`, every selected set of up to 512 output positions has the same exact uniform distribution when the GF seed ranges uniformly over all coefficient vectors.

The theorem does not claim that a fixed deterministic stream contains physical randomness.

---

# 6. Word granularity

The exact production theorem is stated for **64-bit output words**:

> Any selection of at most 512 distinct `uint64_t` output positions is jointly uniform.

Arbitrary subsets of bits inside those same words inherit the corresponding uniformity.

This should not be confused with the stronger claim

> “Any arbitrary 32768 bit positions anywhere in the whole stream are jointly independent.”

For example, taking one bit from each of more than 512 different output words lies outside the stated theorem even though fewer than 32768 total bits may have been selected.

---

# 7. Seed modes

## 7.1 Full theoretical mode

```cpp
uint8_t chacha_seed[32];
uint8_t gf_seed[8192];

ChaCha20GF512FFT rng(chacha_seed, gf_seed);
```

The explicit seed representation is

```text
256 ChaCha bits + 65536 GF bits
= 65792 bits
= 8224 bytes.
```

For the exact probabilistic statement:

- the 8192-byte GF seed must be uniformly sampled from the complete `2^65536` coefficient space;
- if the ChaCha key is also treated probabilistically, it should be sampled independently of the GF seed;
- `stream_id` is a public domain-separation parameter, not hidden seed entropy.

For every fixed ChaCha stream, the GF counting theorem already holds.

---

## 7.2 64-bit convenience mode

```cpp
ChaCha20GF512FFT rng(0x123456789ABCDEF0ULL);
```

This constructor expands one 64-bit value with SplitMix64 to initialize both components.

It exists for:

- deterministic tests;
- examples;
- benchmarks;
- reproducible experiments.

It creates at most `2^64` complete initialized generators, so the full information-theoretic theorem over the explicit GF seed space does **not** apply to this convenience mode.

No deterministic expansion can create 65536 bits of information-theoretic seed entropy from 64 input bits.

---

## 7.3 Practical OS-backed initialization

`ChaCha20GF512Seed.h` provides practical Windows/Linux initialization.

It automatically uses

```cpp
ChaCha20GF512FFT::FULL_SEED_BYTES
```

so it currently reads exactly 8224 bytes.

Windows:

```text
BCryptGenRandom(..., BCRYPT_USE_SYSTEM_PREFERRED_RNG)
```

Linux:

```text
getrandom()
```

Example:

```cpp
#include "ChaCha20GF512Seed.h"

auto seed = chacha20gf512_seed::make_os_full_seed();
ChaCha20GF512FFT rng(seed.data());
```

The Linux helper loops until the complete buffer is filled, because `getrandom()` may legally return a partial result for large requests.

OS-backed CSPRNG output is the recommended practical initialization method.

It should not be confused with a claim that calling the OS API injected 65792 bits of fresh physical entropy. The exact theorem is a statement about uniform sampling of the explicit coefficient space; OS-backed seeding is the practical computational approximation used by normal software.

The supplied `main.cpp` also contains an optional x86 `RDRAND` seeding demonstration. `RDRAND` is not part of the core generator or its formal seed model.

---

# 8. Fast evaluation with an additive FFT

Direct Horner evaluation of a degree-511 polynomial would require

```text
511 GF multiplications per field evaluation.
```

The production implementation instead evaluates 512 GF points at once with a 9-stage additive FFT. Those 512 evaluations supply **1024 output words** because both 64-bit halves are used.

For an aligned GF evaluation block

```text
B = 512*k,
```

the field points are

```text
B, B+1, ..., B+511.
```

The corresponding output-word block is

```text
2B, 2B+1, ..., 2B+1023.
```

Because the lower nine bits of `B` are zero,

```text
B + j = B XOR j
```

for `0 <= j < 512`. Addition in `GF(2^128)` is XOR, so each aligned evaluation block is an affine 9-dimensional binary subspace.

The 512 monomial coefficients are converted once during initialization into a normalized subspace/novel basis. Steady-state evaluation then uses the additive FFT.

The field-operation count per 512 evaluations remains

```text
9 FFT stages * 256 butterfly multiplications = 2304
9 affine-base evaluations                    =   90
----------------------------------------------------
Total                                         = 2394
```

but those 2394 multiplications now produce 1024 output words, or about

```text
2.34 GF multiplications per output word.
```

On x86/x86-64 the fast path uses `PCLMULQDQ`. GCC/Clang now dispatch once per FFT block into a PCLMUL-targeted butterfly routine, so the hot loop no longer crosses the target-attribute boundary for every multiplication. The portable shift/XOR reference path remains available and bit-identical.

The FFT is only an evaluation optimization. It computes the same polynomial values as direct Horner evaluation.

---

# 9. Random access and parallelism

Both components are position indexed.

The main API accepts either the full position type

```cpp
ChaCha20GF512FFT::Position p{low64, high64};
rng.seek(p);
```

or the backward-compatible 64-bit overload. A seek does not need any previous outputs, and seeking within the currently cached 1024-word GF block no longer discards that FFT result.

`ChaCha20GF512FFT` also provides

```cpp
rng.generate(dst, count);
```

for bulk serial generation. The ChaCha component writes full blocks directly into the destination buffer, and the GF layer XORs cached low/high pairs in bulk.

`ChaCha20GF512Parallel.h` provides parallel bulk generation over the **same logical stream**. It does not create multiple RNG streams. Instead:

1. the requested position interval is split into disjoint ranges;
2. worker boundaries are aligned to 1024-output-word / 512-evaluation FFT blocks where possible;
3. each worker copies the prototype generator and seeks to its own start position;
4. all worker results are written into their original places in the output buffer;
5. worker exceptions are captured and rethrown in the calling thread rather than escaping through `std::thread`.

Therefore

```text
parallel output == serial output
```

bit-for-bit, including across the `2^64` boundary of the 128-bit output position.

---

# 10. Performance

The main performance changes in this revision are:

- one 512-point GF FFT now supplies 1024 output words by using both halves;
- GCC/Clang PCLMUL butterflies stay inside one target-specific hot loop;
- ChaCha buffering uses native `uint64_t` words rather than bytewise store/load reconstruction;
- a bulk `generate()` API removes most per-word bookkeeping in bulk/parallel use;
- seeking within the same GF block preserves the FFT cache.

A representative run in the current test environment (AMD EPYC 9V74, GCC 14.2, `-O3`) with 30 million words produced approximately:

| Generator | M `uint64_t`/s | MiB/s |
|---|---:|---:|
| `ChaCha20Counter128` | 56.8 | 433 |
| `ChaCha20GF512 FFT single` (`next_int`) | 26.5 | 203 |
| `ChaCha20GF512 FFT bulk` | 30.3 | 231 |
| `ChaCha20GF512 FFT parallel x4` | 98.5 | 751 |

A separate tight single-thread microbenchmark on the same machine measured the old revision at about 17 M words/s and the current revision at about 31 M words/s. Treat these as implementation measurements, not universal performance claims: CPU frequency, virtualization, compiler, PCLMUL throughput and thread scheduling matter substantially.

---

# 11. Memory footprint

The current implementation stores 128-bit field elements.

Major 512-element arrays are therefore:

```text
monomial coefficient copy   512 * 16 = 8192 bytes
novel-basis coefficient set 512 * 16 = 8192 bytes
FFT result cache            512 * 16 = 8192 bytes
```

plus:

- subspace-polynomial tables;
- normalization constants;
- ChaCha state and block buffer;
- position/cache bookkeeping.

A current build reports

```text
sizeof(ChaCha20GF512FFT) = 27664 bytes.
```

The retained monomial copy is used by the reference Horner path and also keeps one stable header-only class layout in all translation units.

These sizes are implementation details, not part of the mathematical theorem.

---

# 12. Validation

The supplied `main.cpp` checks the current production implementation in several independent ways. A representative successful run includes:

```text
PCLMUL runtime support:    yes
sizeof(ChaCha20GF512FFT):  27664 bytes
Full seed size:            8224 bytes
GF seed size:              8192 bytes
ChaCha128 regression vec:  OK
Hybrid regression vector: OK
GF(2^128) fast/ref:        OK
GF128 FFT vs Horner:       OK
Hybrid 128-bit seek:       OK
Bulk vs next_int:          OK
Parallel identity x4:      OK
```

The important implementation checks are:

### ChaCha regression

A fixed vector makes unintended changes to the position-indexed ChaCha component visible.

### Full hybrid regression

A second known-answer vector covers the complete `ChaCha XOR GF` stream. This catches accidental stream changes that a GF-vs-Horner test alone cannot detect.

### GF multiplication

The portable `GF(2^128)` multiplication is checked against the PCLMUL implementation.

### FFT vs. Horner

The additive FFT is compared bit-for-bit with direct polynomial evaluation, including positions with a nonzero high 64-bit limb.

### 128-bit seek, bulk and serial/parallel identity

Random access is checked against sequential generation across the low-64-bit carry. The bulk API is compared against repeated `next_int()` calls at ordinary positions, 1024-word GF boundaries, across the low-64-bit carry and exactly at the final position `2^128-1`. Parallel generation is checked against the same serial logical stream.

### Toy-model theorem test

`toy_kwise_test.h` exhaustively checks scaled-down versions of the construction:

```text
(A) k-wise independence
(B) sharpness of the unprojected polynomial construction at k+1
(C) invariance under XOR with an arbitrary fixed position-dependent mask
(D) negative control with a reducible polynomial
(E) projection from a larger field element to a smaller output word
(F) splitting each field value into low/high adjacent output halves
```

The new `(F)` check mirrors the production change with `GF(2^4)`: each 4-bit field evaluation is split into two 2-bit output words, and mixed selections (including both halves of one evaluation) are verified to remain exactly uniform.

The demo uses the quick toy mode by default. Passing `full` as the third command-line argument after count and thread count runs the complete GF(2^4) checks plus the sampled GF(2^8) checks. The GF(2^8) sharpness test now stores only the `2^24` generated 32-bit tuples (64 MiB) rather than allocating a 512-MiB bitmap for the full `2^32` tuple universe.

A full run of the revised tests passes.

---

# 13. Portability and build

The project is C++17 and is intended for:

- Visual Studio 2022 / MSVC;
- GCC;
- Clang.

On x86/x86-64 it detects `PCLMULQDQ` at runtime.

GCC/Clang use a function-specific target attribute together with

```cpp
__builtin_cpu_supports("pclmul")
```

so the complete program does not need to be compiled globally with `-mpclmul`.

The MSVC path uses the corresponding x86/x64 intrinsics.

## Direct GCC build

```bash
g++ -std=c++17 -O3 -pthread main.cpp -o chacha20gf512
./chacha20gf512 1000000 4
```

To force the portable GF path:

```bash
g++ -std=c++17 -O3 -pthread -DCHACHA20GF512_FORCE_PORTABLE \
    main.cpp -o chacha20gf512-portable
```

The portable path must produce the same stream; it is only slower.

## CMake

If the supplied `CMakeLists.txt` is used:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

On Visual Studio generators, `--config Release` selects the optimized configuration.

---

# 14. API examples

## Full explicit seed

```cpp
#include "ChaCha20GF512FFT.h"

uint8_t chacha_seed[ChaCha20GF512FFT::CHACHA_SEED_BYTES] = {};
uint8_t gf_seed[ChaCha20GF512FFT::GF_SEED_BYTES] = {};

ChaCha20GF512FFT rng(chacha_seed, gf_seed);

uint64_t a = rng.next_int();
uint64_t b = rng.next_int();

std::vector<uint64_t> block(4096);
rng.generate(block.data(), block.size());
```

## 128-bit random access

```cpp
ChaCha20GF512FFT::Position p{
    0x0123456789abcdefULL,   // low 64 bits
    0x0000000000000001ULL    // high 64 bits
};

rng.seek(p);
uint64_t x = rng.next_int();
```

## Backward-compatible 64-bit seek

```cpp
rng.seek(1'000'000ULL);
```

## Practical OS-backed seed

```cpp
#include "ChaCha20GF512Seed.h"

auto seed = chacha20gf512_seed::make_os_full_seed();
ChaCha20GF512FFT rng(seed.data());
```

## Parallel generation of the same stream

```cpp
#include "ChaCha20GF512Parallel.h"

ChaCha20GF512Parallel prng(seed.data());

std::vector<uint64_t> values(1'000'000);
prng.fill_parallel(values.data(), values.size(), 4);
```

The parallel wrapper also accepts the 128-bit `Position` type for `seek()` and `fill_parallel_at()`.

---

# 15. What is *not* claimed

ChaCha20GF512 does **not** claim:

- to be empirically “more random than ChaCha20”;
- to replace established cryptographic standards;
- to be the fastest or smallest PRNG;
- unlimited output;
- that the 64-bit output values have a proven minimal period of `2^128`;
- that output words never repeat;
- exact independence beyond the proven 512-word order;
- that 513-wise independence of the split-half 64-bit stream has been proved false;
- independence of arbitrary bit sets spread across more than 512 output words;
- that the GF component increases ChaCha20's cryptographic key strength;
- that an OS call returning 8224 bytes injected 65792 bits of fresh physical entropy;
- that the 64-bit convenience constructor obtains the full-seed theorem;
- backtracking resistance after compromise of a stateful DRBG state;
- cryptographic audit or suitability as a drop-in system CSPRNG;
- novelty of ChaCha20, HChaCha20, polynomial `k`-wise independence, additive FFTs, GHASH-style field arithmetic or XOR-combiner reasoning;
- universal mathematical optimality across every conceivable RNG property.

The phrase

> **“Maximizing What Can Be Proven About a PRNG”**

describes the design motivation, not a theorem of universal optimality.

---

# 16. FAQ / likely objections

## “Why use GF(2^128) if the public output is 64-bit words?”

Because the polynomial construction needs a large field for both coefficients and interpolation points, and `GF(2^128)` is especially convenient on modern CPUs because carry-less multiplication is directly accelerated.

The current mapping also uses the field width completely: each uniform 128-bit evaluation becomes two adjacent 64-bit output words. A `2^128`-word stream therefore consumes `2^127` distinct GF evaluation points, still far below the `2^128` available points.

---

## “Is 2^128 the period?”

Not in the traditional recurrence-generator sense.

The implementation defines a non-wrapping position domain of exactly `2^128` words. Every position is unique, and generation stops rather than reusing position zero.

No claim is made that the sequence of 64-bit values has minimal period `2^128`.

The practical requirement motivating the change was simpler:

> **period exhaustion or counter wrap should be impossible in any physically meaningful use.**

A `2^128`-word non-wrapping stream meets that requirement by an enormous margin.

---

## “Why not add an even larger-period LFSR?”

A separate maximal-period component could make a larger formal cycle, but the exact polynomial guarantee would still be tied to its finite set of distinct field positions.

The current design instead makes the exact theorem valid over the complete supported stream and then stops cleanly at the end of that domain.

That is simpler to state and prove.

---

## “Does splitting 128 bits into two 64-bit words waste randomness?”

No field bits are discarded now.

For each evaluation `P(j)`, output position `2j` uses the low 64 bits and `2j+1` uses the high 64 bits. Together those two words contain the complete 128-bit field value.

If only one half is selected in a theorem query, it is a uniform 64-bit projection. If both halves are selected, they are jointly the original uniform 128-bit value.

---

## “Why does the GF seed double?”

A degree-511 polynomial still has 512 coefficients, but each coefficient is now a 128-bit field element:

```text
512 * 128 = 65536 bits.
```

The extra seed space pays for the larger field and therefore the larger distinct-position domain.

---

## “Is the new seed size minimal?”

Not for the statement “arbitrary 512 output words are independent.” That theorem has a lower bound of 32768 GF seed bits, while this construction uses 65536.

However, the current split-half mapping no longer throws away the extra dimension. For 512 distinct evaluation points with both halves exposed (1024 output words), the 65536 output bits are exactly uniform and use the complete 65536-bit GF seed space one-to-one.

That structured fact must not be confused with arbitrary 1024-wise independence.

---

## “Could the 64-bit output actually be more than 512-wise independent?”

Possibly.

The stream now exposes low and high coordinates at adjacent positions, but a set of more than 512 words may still choose only one coordinate from each of more than 512 distinct field evaluations. Those projections can hide relations that exist between the full 128-bit evaluations.

This project therefore proves and tests the conservative statement that **512-wise independence definitely holds**, while seed counting gives an absolute upper bound of 1024. The exact maximum order in between remains open here.

---

## “Why use HChaCha20?”

A full 128-bit ChaCha block index needs all four state words 12..15.

HChaCha20 moves stream separation into the derived key instead:

```text
subkey = HChaCha20(key, stream_id || 0)
```

so words 12..15 remain completely available for position addressing.

---

## “Why not derive the 8192-byte GF seed from the ChaCha key?”

Because deterministic expansion of a 256-bit key reaches at most `2^256` GF coefficient vectors.

The exact information-theoretic theorem stated here is a counting theorem over the complete explicit GF coefficient space.

A KDF/CSPRNG expansion from a short seed may be computationally excellent, but it does not create the missing information-theoretic seed space.

---

## “Isn't the GF part algebraically weak?”

Yes, deliberately.

Bounded independence and cryptographic pseudorandomness are different properties.

The GF polynomial is included because its finite-order distribution property is exact and transparent. ChaCha is included to provide the computationally pseudorandom mask.

---

## “Why exactly 512?”

`k = 512` remains a useful practical/theoretical point:

- 512 output words cover 32768 output bits;
- the polynomial degree is still manageable;
- 512 evaluation points map naturally to a 9-stage additive FFT;
- an 8192-byte GF coefficient set is still small in ordinary software;
- the parallel implementation can split naturally at 1024-output-word boundaries, each backed by one 512-point FFT.

Other values of `k` are mathematically possible.

---

# 17. Repository layout

The current repository consists of:

```text
README.md
CMakeLists.txt
ChaCha20GF512FFT.h
ChaCha20GF512Parallel.h
ChaCha20GF512Seed.h
toy_kwise_test.h
main.cpp
```

Roles:

- `ChaCha20GF512FFT.h` — generator, GF arithmetic, additive FFT, ChaCha/HChaCha position logic;
- `ChaCha20GF512Parallel.h` — bulk parallel generation of the same logical stream;
- `ChaCha20GF512Seed.h` — practical Windows/Linux OS-backed seed acquisition;
- `toy_kwise_test.h` — exhaustive small-field verification of the mathematical construction;
- `main.cpp` — validation, seed demonstrations and benchmarks;
- `CMakeLists.txt` — optional cross-platform build.

The core generator itself does not depend on the OS random APIs or on the parallel wrapper.

---

# 18. Future work

Possible directions include:

- VPCLMUL / wider SIMD acceleration of GF(2^128) FFT work;
- ARM `PMULL` acceleration;
- reducing the cost of 128-bit field multiplication;
- more benchmark data across CPUs and compilers;
- PractRand, TestU01/BigCrush and SmokeRand as implementation sanity checks rather than proofs;
- parameterizing the independence order `k`;
- studying the **exact independence order of the split-half 64-bit family above 512**;
- alternative linear output mappings from each 128-bit field evaluation;
- independent cryptographic-combiner variants if robustness against failure of one cryptographic primitive becomes a separate design goal.

These are possible extensions, not claims of the current implementation.

---

# 19. References

## ChaCha20

Daniel J. Bernstein, **“ChaCha, a variant of Salsa20”**, 2008.  
https://cr.yp.to/chacha/chacha-20080120.pdf

Y. Nir and A. Langley, **RFC 8439: ChaCha20 and Poly1305 for IETF Protocols**.  
RFC 8439 documents the IETF layout and explicitly notes that original ChaCha used a 64-bit nonce and 64-bit block counter.  
https://www.rfc-editor.org/rfc/rfc8439.html

## HChaCha20 / XChaCha

S. Arciszewski, **XChaCha: eXtended-nonce ChaCha and AEAD_XChaCha20_Poly1305**, Internet-Draft.  
Contains the HChaCha20 construction and output-word selection `0,1,2,3,12,13,14,15`.  
https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-xchacha

## k-wise independence

James Aspnes, **“k-wise Independence”**.  
Describes the standard random degree-`k-1` polynomial construction over a finite field.  
https://www.cs.yale.edu/homes/aspnes/pinewiki/KwiseIndependence.html

## GF(2^128)

NIST SP 800-38D, **Recommendation for Block Cipher Modes of Operation: Galois/Counter Mode (GCM) and GMAC**.  
Defines the 128-bit binary field used by GHASH.  
https://csrc.nist.gov/pubs/sp/800/38/d/final

## Additive FFT

Shuhong Gao and Todd Mateer, **“Additive Fast Fourier Transforms over Finite Fields”**, IEEE Transactions on Information Theory 56(12), 2010.  
https://www.math.clemson.edu/~sgao/papers/GM10.pdf

Sian-Jheng Lin, Wei-Ho Chung, and Yunghsiang S. Han, **“Novel Polynomial Basis and Its Application to Reed-Solomon Erasure Codes”**.  
https://arxiv.org/abs/1404.3458

## Carry-less multiplication

Intel, **“Enabling High-Performance Galois-Counter-Mode on Intel Architecture Processors”**.  
Discusses 128-bit binary-field multiplication with `PCLMULQDQ`, including classical and Karatsuba variants.  
https://www.intel.com/content/dam/www/public/us/en/documents/software-support/enabling-high-performance-gcm.pdf

## OS-backed seed acquisition

Linux `getrandom(2)` manual page.  
https://man7.org/linux/man-pages/man2/getrandom.2.html

Microsoft, **BCryptGenRandom**.  
https://learn.microsoft.com/en-us/windows/win32/api/bcrypt/nf-bcrypt-bcryptgenrandom

---

# Status

The current implementation is a research-quality experimental prototype.

The 128-bit revision has independent checks for:

- portable vs. PCLMUL GF(2^128) arithmetic;
- additive FFT vs. direct Horner evaluation;
- positions above `2^64`;
- carry across the 64-bit position boundary;
- serial vs. parallel identity;
- exhaustive toy-model `k`-wise independence;
- XOR masking;
- a negative reducible-field control;
- and the low/high split-field output mapping.

The design goal is to make every major claim separable:

- what is exact;
- what is computational;
- what is an implementation optimization;
- what is merely a practical initialization choice;
- and what is not claimed at all.

Feedback, counterexamples, theoretical corrections, implementation review and comparisons with related constructions are welcome.
