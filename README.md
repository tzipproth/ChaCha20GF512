# ChaCha20GF512
## Maximizing What Can Be Proven About a PRNG

> **Once a generator already passes every practical statistical test, what theoretical properties can still be improved?**

ChaCha20GF512 is an experimental pseudorandom number generator built around that question.

The goal is **not** to claim that its output is empirically “more random than ChaCha20.” Modern high-quality generators can reach a point where practical statistical test suites no longer provide a useful ranking of their randomness. ChaCha20GF512 therefore takes a different approach: it combines two complementary notions of pseudorandomness and asks how much can be said about the result **exactly**, rather than only empirically.

The construction combines:

1. **ChaCha20**, used in its original 256-bit-key / 64-bit-counter / 64-bit-stream-ID layout, as the computational pseudorandom component.
2. A random degree-511 polynomial over the **Galois field GF(2^64)**, supplying exact **512-wise independence** for 64-bit output words.

**GF** stands for **Galois field**, another name for a finite field. `GF(2^64)` contains exactly `2^64` elements, and in this implementation each field element is represented by one 64-bit word. Field addition is bitwise XOR; multiplication is carry-less polynomial multiplication followed by reduction modulo a fixed irreducible degree-64 polynomial. This field structure is what makes polynomial interpolation work and, in turn, gives the exact 512-wise independence result used by ChaCha20GF512.

The final output is simply

```text
R(i) = ChaCha20(K, stream_id, i) XOR P(i)
```

where

```text
P(x) = a0 + a1*x + ... + a511*x^511    over GF(2^64).
```

The central claim is therefore not

> “ChaCha20GF512 is more random than ChaCha20.”

but rather

> **ChaCha20GF512 makes strictly stronger information-theoretic statements about its output family while retaining the ordinary computational pseudorandomness argument of ChaCha20.**

This is a research/experimental PRNG, not a cryptographic standard and not a claim of a universally “optimal” RNG.

---

## Why this is a genuine theoretical addition

With a fixed/public `stream_id`, ordinary ChaCha20 is indexed by a 256-bit key. Therefore there are at most

```text
2^256
```

distinct keyed streams.

Exact 5-wise independence for 64-bit output words would require five selected outputs to cover all

```text
2^(5*64) = 2^320
```

joint values with equal probability.

A family containing only `2^256` streams cannot do that. So, under this model, **ChaCha20 alone cannot even be exactly 5-wise independent for 64-bit words**. This is an information-theoretic counting argument; it says nothing negative about ChaCha20's cryptographic strength.

ChaCha20GF512 adds a separate 32768-bit GF seed space and reaches exact 512-wise independence:

```text
512 * 64 = 32768 bits.
```

Thus the hybrid possesses an exact distribution property that a 256-bit-key ChaCha20 family **provably cannot possess**, while the XOR construction does not sacrifice ChaCha20's computational pseudorandomness under the assumptions stated below.

This is the main motivation for the project.

---

# 1. Construction

## 1.1 ChaCha20 component

The implementation uses the original ChaCha layout:

- 256-bit key;
- 64-bit block counter;
- 64-bit stream ID / nonce field;
- 20 rounds.

RFC 8439 standardizes the later IETF 32-bit-counter / 96-bit-nonce layout, but explicitly notes that the original ChaCha construction used a **64-bit nonce and 64-bit block counter**.

ChaCha20GF512 uses ChaCha20 as a counter-indexed stream primitive, not as a state-updating DRBG. There is no periodic rekeying. That is intentional: backtracking resistance after internal-state compromise is a separate DRBG property, whereas this project is concerned with pseudorandomness, exact distribution guarantees, random access, and parallel evaluation.

For a fixed key and stream ID, each ChaCha block can be addressed directly by its counter.

## 1.2 Galois field GF(2^64) component

The second component is a degree-at-most-511 polynomial

```text
P(x) = a0 + a1*x + ... + a511*x^511
```

over the field `GF(2^64)`.

There are 512 coefficients, each one 64-bit field element:

```text
512 * 64 bits = 32768 bits = 4096 bytes.
```

The implementation represents `GF(2^64)` modulo the primitive polynomial

```text
x^64 + x^4 + x^3 + x + 1
```

with reduction constant `0x1B`.

The mathematical GF output at word position `i` is

```text
G(i) = P(i),
```

where the 64-bit integer `i` is interpreted bijectively as a field element.

## 1.3 Combination

The final output word is

```text
R(i) = C(i) XOR G(i),
```

where `C(i)` is the corresponding 64-bit ChaCha20 output.

XOR is useful here because XOR with a fixed vector is a bijection. It therefore preserves exact uniformity of the GF component, and it also gives a simple reduction argument for preservation of ChaCha20's computational pseudorandomness when the two seed components are independent.

---

# 2. What is proven exactly

## 2.1 Exact 512-wise independence

Let

```text
P(x) = a0 + a1*x + ... + a511*x^511
```

where all 512 coefficients are independent and uniformly distributed in `GF(2^64)`.

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
GF(2^64)^t.
```

The reason is the Vandermonde/interpolation argument: prescribing values at `t` distinct points imposes `t` independent linear constraints on the 512 coefficients. For `t = 512`, there is exactly one degree-at-most-511 polynomial producing any chosen 512-tuple.

Therefore, for **any 512 distinct output-word positions**, all

```text
2^(512*64) = 2^32768
```

possible 512-tuples occur exactly once over the complete 4096-byte GF seed space.

For `t < 512`, every possible `t`-tuple occurs exactly

```text
2^(32768 - 64*t)
```

times over that seed space.

This is an exact combinatorial statement, not a statistical-test result.

## 2.2 Arbitrary positions, not only consecutive outputs

The selected positions need not be adjacent. For example:

```text
0
17
4711
123456789
...
```

are just as valid as 512 consecutive positions.

The only requirement is that the selected 64-bit positions correspond to distinct field points. Since `uint64_t -> GF(2^64)` is bijective, this holds for distinct positions in

```text
0 <= i < 2^64.
```

That domain contains `2^64` output words, or 128 EiB of output data.

## 2.3 XOR with a fixed ChaCha stream preserves the exact guarantee

Fix any ChaCha key and any set of `t <= 512` output positions. At those positions ChaCha contributes some fixed vector

```text
C = (C1, ..., Ct).
```

The GF component contributes an exactly uniform vector `G`.

The hybrid produces

```text
R = C XOR G.
```

Because XOR with a fixed vector is a bijection, `R` is exactly uniform too.

So the 512-wise guarantee does **not** depend on ChaCha20 being cryptographically secure. Even a completely known fixed mask cannot destroy it.

If the ChaCha key is also random, the same probabilistic guarantee holds when the ChaCha key and GF seed are sampled independently.

## 2.4 Information-theoretically minimal GF seed size

Exact 512-wise independence of 64-bit words requires all

```text
2^32768
```

possible joint outcomes of 512 selected words to be representable with equal probability.

Therefore any exact construction needs at least

```text
32768 seed bits.
```

The GF component uses exactly 32768 bits.

Thus, **with respect to seed-space size, the GF construction reaches the information-theoretic lower bound for exact 512-wise independence of 64-bit words.**

---

# 3. Where the exact guarantee stops

The GF component is deliberately simple and algebraic. It is linear in its 512 field coefficients.

For any fixed set of **513 distinct positions**, the map

```text
(a0, ..., a511) -> (P(x1), ..., P(x513))
```

maps a 512-dimensional vector space into a 513-dimensional one. Its image cannot fill the entire output space.

Consequently, for every chosen set of 513 positions there exists a non-trivial GF-linear relation

```text
lambda1*P(x1) + ... + lambda513*P(x513) = 0
```

that holds for **every** GF seed; the coefficients depend only on the selected positions.

This is exactly why the GF component alone is not a general-purpose cryptographic PRG. Once enough unmasked evaluations are known, the polynomial can be reconstructed by interpolation.

The intended division of labor is therefore:

| Component | What it provides | Where it stops |
|---|---|---|
| GF polynomial | Exact information-theoretic independence through 512 output words | Explicit algebraic structure exists beyond that order |
| ChaCha20 | Computational pseudorandomness with no known efficient distinguisher | No comparable exact high-order independence follows from a 256-bit key space |
| ChaCha20GF512 | Both guarantees coexist | Larger state and lower throughput than ChaCha20 alone |

This boundary is part of the design, not something hidden by it.

---

# 4. What relies on the ChaCha20 assumption

The exact 512-wise result above is unconditional once the GF seed is sampled as specified.

The claim beyond that finite exact order is different: it relies on ChaCha20 being computationally pseudorandom.

Let `C` be a ChaCha20 stream and `G` an independently generated GF stream. Suppose an efficient distinguisher could distinguish

```text
C XOR G
```

from a uniform random stream.

Then the same distinguisher could be used against ChaCha20: given a challenge stream `X` that is either ChaCha20 or uniform random, independently sample `G` and give

```text
X XOR G
```

to the distinguisher.

If `X` is uniform random, `X XOR G` is still uniform random. If `X` is ChaCha20, the result is the hybrid construction. Therefore a distinguisher for the hybrid would imply a distinguisher for ChaCha20.

This is a reduction argument, not a proof that ChaCha20 itself is secure.

## Graceful degradation of guarantees

The two guarantees can fail separately:

- if a future cryptanalytic breakthrough invalidated ChaCha20's computational pseudorandomness assumption, the **exact 512-wise GF guarantee would remain**;
- if one ignores the GF theorem beyond 512 words, the ChaCha component still supplies the usual computational pseudorandomness argument.

This is one reason for combining two structurally different components.

It does **not** mean the hybrid would remain cryptographically secure after a hypothetical break of ChaCha20; it means only that the information-theoretic 512-wise statement survives.

---

# 5. The guarantee is over a generator family

A PRNG instance with a fixed seed is deterministic. After the seed is fixed, there is no probability distribution left inside that one stream.

The phrase

> “ChaCha20GF512 is 512-wise independent”

is therefore shorthand for a statement about the **family of streams obtained by sampling the GF coefficients uniformly**.

Equivalently, it is a combinatorial statement over the complete GF seed space.

This distinction matters. ChaCha20GF512 does not somehow turn a deterministic sequence into physical randomness; it constructs a deterministic family with unusually strong exact distribution properties over its initialization space.

---

# 6. Word granularity

The exact theorem is stated for **64-bit output words**:

> Any selection of at most 512 distinct `uint64_t` output positions is jointly uniform.

Since a uniform vector of 512 words contains 32768 uniform bits, arbitrary subsets of bits *inside those same 512 words* inherit the corresponding uniformity.

But this should not be confused with the stronger claim

> “Any arbitrary 32768 bit positions anywhere in the whole stream are jointly independent.”

For example, 513 single bits taken from 513 different output words are already outside the proven 512-word guarantee.

---

# 7. Seed modes

## 7.1 Full theoretical mode

```cpp
uint8_t chacha_seed[32];
uint8_t gf_seed[4096];

ChaCha20GF512FFT rng(chacha_seed, gf_seed);
```

The complete explicit seed representation is

```text
256 ChaCha bits + 32768 GF bits = 33024 bits = 4128 bytes.
```

For the full probabilistic statement:

- the 4096-byte GF seed must be uniformly sampled from the full `2^32768` coefficient space;
- if the ChaCha key is also random, it should be sampled independently of the GF seed;
- `stream_id` is treated as a fixed/public domain-separation parameter, not as hidden seed entropy.

There are two related statements:

**Combinatorial:** for every fixed ChaCha stream, every selected `t <= 512` hybrid output tuple has exactly the same number of GF-seed preimages.

**Probabilistic:** if the GF seed is uniform, those selected outputs are exactly jointly uniform random variables.

## 7.2 64-bit convenience mode

```cpp
ChaCha20GF512FFT rng(0x123456789ABCDEF0ULL);
```

This constructor expands one 64-bit value with SplitMix64 to initialize both components.

It exists only for reproducibility, examples, and benchmarks.

It creates at most

```text
2^64
```

complete initialized generators. Therefore:

- the full 512-wise information-theoretic guarantee does **not** apply;
- the ChaCha key has at most 64 bits of originating seed entropy, not 256;
- the independent-seed reduction argument does not apply formally because both components come from the same short seed.

No deterministic expansion can create additional seed entropy or additional reachable initial states.

### Do I really need 4096 bytes of independent seed material?

For the **exact probabilistic 512-wise statement**, yes: the GF coefficient vector must be sampled uniformly from its full `2^32768` space. This is a mathematical initialization model, not a claim that applications routinely collect 32768 bits of fresh physical entropy.

If the 4096-byte GF seed is instead expanded deterministically from a shorter seed using a CSPRNG or KDF, the resulting stream may still be computationally excellent, but the exact information-theoretic 512-wise claim over initialization no longer follows. The distinction is intentional.

## 7.3 Practical OS-backed initialization

`ChaCha20GF512Seed.h` is an optional convenience header for normal applications. It obtains the complete 4128-byte explicit seed from the operating system:

- Windows: `BCryptGenRandom` with `BCRYPT_USE_SYSTEM_PREFERRED_RNG`;
- Linux: `getrandom()`.

Example:

```cpp
#include "ChaCha20GF512Seed.h"

auto seed = chacha20gf512_seed::make_os_full_seed();
ChaCha20GF512FFT rng(seed.data());
```

This is the recommended **practical** initialization mode, but it is deliberately not advertised as 33024 bits of fresh physical entropy. Operating-system random APIs are CSPRNG interfaces. The information-theoretic theorem still refers to uniform sampling of the full explicit GF seed space.

The supplied `main.cpp` also contains an optional x86 `RDRAND` seeding demonstration. That code is intentionally kept out of the generator and seed-helper APIs: `RDRAND` is hardware-DRBG output and is shown only as an example, not as part of the ChaCha20GF512 construction or its entropy claim.

---

# 8. Fast evaluation with an additive FFT

The first implementation evaluated the degree-511 polynomial independently at each point with Horner's rule:

```text
511 GF multiplications per output word.
```

That was mathematically simple but expensive.

The current implementation computes **exactly the same polynomial values** with a 9-stage additive FFT.

For every aligned block

```text
B = 512*k,
```

the 512 evaluation points are

```text
{ B, B+1, ..., B+511 }
```

and, because the lower nine bits of `B` are zero,

```text
B + j = B XOR j    for 0 <= j < 512.
```

Addition in `GF(2^64)` is XOR, so these points form the affine subspace

```text
B + span_F2(1, 2, 4, ..., 256).
```

This is precisely the structure exploited by additive FFTs over characteristic-two finite fields.

The 512 monomial coefficients supplied by the seed are converted once during initialization into a normalized subspace/novel polynomial basis. Steady-state evaluation then uses the additive FFT.

Current multiplication count per 512-word GF block:

```text
9 stages * 256 butterfly multiplications = 2304
9 affine-base evaluations                 =   90
------------------------------------------------
Total                                      = 2394
```

or approximately

```text
4.68 GF multiplications per output word
```

instead of 511 with direct Horner evaluation.

The FFT is only an evaluation optimization. **It does not change the mathematical generator.** The supplied test program verifies FFT output bit-for-bit against direct Horner evaluation at fixed and random positions.

---

# 9. Random access and parallelism

Both components are position-indexed:

```text
ChaCha block = ChaCha20(key, block_counter)
GF word      = P(word_position)
```

so previous outputs are not needed to compute a later one.

The API provides

```cpp
rng.seek(position);
```

and independent 512-word ranges can be evaluated independently, making parallel evaluation conceptually straightforward.

The FFT implementation caches one aligned 512-word GF block. A seek to another block recomputes only the target block.

`ChaCha20GF512Parallel.h` adds a bulk parallel convenience wrapper:

```cpp
#include "ChaCha20GF512Parallel.h"

ChaCha20GF512Parallel rng(full_seed);
std::vector<uint64_t> values(1'000'000);
rng.fill_parallel(values.data(), values.size(), 4);
```

This does **not** create or combine multiple random streams. Workers evaluate disjoint position ranges of the same logical stream. Range splitting is aligned to 512-word FFT blocks where possible, so the parallel result is bit-identical to serial generation while avoiding duplicated FFT work at worker boundaries.

A `fill_parallel_at()` overload provides position-indexed bulk generation, and a thread count of zero selects `std::thread::hardware_concurrency()`. Thread creation is intentionally outside the core RNG; the wrapper is aimed at large bulk requests, not individual `next_int()` calls.

---

# 10. Performance

The current implementation was benchmarked on two x86-64 CPUs under both
Windows/MSVC and Linux/GCC. All runs below used one million `uint64_t`
outputs, `PCLMULQDQ` runtime support, and the same deterministic benchmark
seed/stream. The single-threaded and parallel hybrid runs produce the same
logical output stream.

## Summary

| CPU / toolchain | ChaCha20Counter64 | ChaCha20GF512 FFT single | ChaCha20GF512 FFT parallel x4 |
|---|---:|---:|---:|
| Intel Core i7-4790K / Windows MSVC | 41.160 M/s | 18.769 M/s | 68.618 M/s |
| Intel Core i7-4790K / Linux GCC | 76.581 M/s | 35.020 M/s | 77.172 M/s |
| AMD Ryzen 9 5900X / Windows MSVC | 53.550 M/s | 31.630 M/s | 105.258 M/s |
| AMD Ryzen 9 5900X / Linux GCC | 99.554 M/s | 43.846 M/s | 140.903 M/s |

The compiler/toolchain has a substantial effect on this implementation,
especially on the scalar ChaCha20 path. The benchmark numbers should therefore
be read as measurements of these particular builds, not as universal ratios.

## Intel Core i7-4790K (Haswell) — Windows / MSVC

```text
PCLMUL runtime support:    yes
sizeof(ChaCha20GF512FFT):    14016 bytes (test build)
ChaCha 64/64 test vector:  OK
GF arithmetic fast/ref:    OK
FFT vs Horner:             OK
Hybrid seek test:          OK
Parallel identity x4:      OK
OS random seed demo:       OK
RDRAND seed demo:          OK

Generator                              M uint64/s       MiB/s
----------------------------------------------------------------
ChaCha20Counter64                         41.160       314.028
ChaCha20GF512 FFT single                    18.769       143.195
ChaCha20GF512 FFT parallel x4               68.618       523.511
```

The single-threaded hybrid reaches about **45.6%** of scalar ChaCha20
throughput, or is roughly **2.19x** slower.

The 4-thread bulk path is about **3.66x** faster than the hybrid single-thread
run, corresponding to roughly **91%** of ideal 4-thread scaling.

## Intel Core i7-4790K (Haswell) — Linux / GCC

Built with:

```bash
g++ -std=c++17 -O3 -pthread main.cpp -o chacha20gf512
```

Measured result:

```text
PCLMUL runtime support:    yes
sizeof(ChaCha20GF512FFT):    14016 bytes (test build)
ChaCha 64/64 test vector:  OK
GF arithmetic fast/ref:    OK
FFT vs Horner:             OK
Hybrid seek test:          OK
Parallel identity x4:      OK
OS random seed demo:       OK
RDRAND seed demo:          OK

Generator                              M uint64/s       MiB/s
----------------------------------------------------------------
ChaCha20Counter64                         76.581       584.270
ChaCha20GF512 FFT single                    35.020       267.182
ChaCha20GF512 FFT parallel x4               77.172       588.779
```

The single-threaded hybrid reaches about **45.7%** of scalar ChaCha20
throughput. The 4-thread path is about **2.20x** faster than the hybrid
single-thread run.

This result is also a useful reminder that WSL/Linux itself is not inherently
a performance obstacle for this CPU-bound code: the GCC build is substantially
faster than the measured MSVC build on the same physical processor.

## AMD Ryzen 9 5900X (Zen 3) — Windows / MSVC

```text
PCLMUL runtime support:    yes
sizeof(ChaCha20GF512FFT):    14016 bytes (test build)
ChaCha 64/64 test vector:  OK
GF arithmetic fast/ref:    OK
FFT vs Horner:             OK
Hybrid seek test:          OK
Parallel identity x4:      OK
OS random seed demo:       OK
RDRAND seed demo:          OK

Generator                              M uint64/s       MiB/s
----------------------------------------------------------------
ChaCha20Counter64                         53.550       408.550
ChaCha20GF512 FFT single                    31.630       241.317
ChaCha20GF512 FFT parallel x4              105.258       803.052
```

The single-threaded hybrid reaches about **59.1%** of scalar ChaCha20
throughput, or is roughly **1.69x** slower.

The 4-thread bulk path is about **3.33x** faster than the hybrid single-thread
run, corresponding to roughly **83%** of ideal 4-thread scaling.

## AMD Ryzen 9 5900X (Zen 3) — Linux / GCC

Built with the same GCC command shown above.

```text
PCLMUL runtime support:    yes
sizeof(ChaCha20GF512FFT):    14016 bytes (test build)
ChaCha 64/64 test vector:  OK
GF arithmetic fast/ref:    OK
FFT vs Horner:             OK
Hybrid seek test:          OK
Parallel identity x4:      OK
OS random seed demo:       OK
RDRAND seed demo:          OK

Generator                              M uint64/s       MiB/s
----------------------------------------------------------------
ChaCha20Counter64                         99.554       759.538
ChaCha20GF512 FFT single                    43.846       334.516
ChaCha20GF512 FFT parallel x4              140.903      1075.005
```

The single-threaded hybrid reaches about **44.0%** of scalar ChaCha20
throughput, or is roughly **2.27x** slower.

The 4-thread bulk path is about **3.21x** faster than the hybrid single-thread
run, corresponding to roughly **80%** of ideal 4-thread scaling. Its measured
throughput is about **1.05 GiB/s**.

## Effect of the additive FFT

The original direct-Horner version of exactly the same GF512 construction
produced only about

```text
1.433 M uint64/s
```

on the i7-4790K Windows benchmark.

The final FFT version on that same Windows/MSVC system reaches

```text
18.769 M uint64/s single-threaded
68.618 M uint64/s with four threads
```

without changing the generated sequence or the mathematical guarantees.

Thus, on that machine, the additive FFT improves the complete hybrid by about
**13.1x** single-threaded, while the 4-thread implementation is about
**47.9x** faster than the original Horner prototype.

## Interpreting the benchmark

These numbers are not intended as universal performance claims. They mainly
show that:

1. exact 512-wise independence is no longer prohibitively expensive once the
   polynomial is evaluated with an additive FFT;
2. both components are position-addressable, so bulk generation parallelizes
   naturally;
3. the parallel implementation computes **the same logical stream** as serial
   generation rather than combining independent per-thread RNG streams;
4. compiler code generation matters substantially for this workload.

The benchmark `sink` value is identical for serial and parallel generation.
That is a useful implementation sanity check: optimization and parallelization
do not change the generated sequence.

CPU architecture, compiler version, optimization flags, clocking, PCLMUL
implementation, thread scheduling, and future SIMD versions can all change
these ratios.

---

# 11. Memory footprint

The mathematical GF seed/state is exactly

```text
4096 bytes.
```

The FFT implementation additionally caches one 512-word evaluated block:

```text
4096 bytes.
```

Together with the transformed coefficients, ChaCha state, FFT plan data, and bookkeeping, the current 64-bit production object occupies approximately

```text
9920 bytes
```

in the tested GCC build.

The supplied test build defines `CHACHA20GF512_ENABLE_TEST_API` and intentionally retains an additional 4096-byte copy of the original monomial coefficients so the FFT can be compared against direct Horner evaluation; that test object occupies 14016 bytes in the current build.

These sizes are implementation details, not part of the mathematical construction.

---

# 12. Validation

The supplied benchmark/test program currently checks:

- a known zero-key / zero-nonce ChaCha20 64/64 test vector;
- portable GF multiplication against the PCLMUL implementation;
- basic field identities such as distributivity and multiplication by one;
- additive FFT results against direct Horner evaluation for multiple unrelated coefficient sets;
- fixed boundary positions and random 64-bit positions;
- random-access `seek()` against sequential generation;
- serial versus parallel bulk generation at an unaligned starting position.

The FFT and Horner paths produce the same hybrid stream bit-for-bit, and the parallel wrapper produces the same logical stream as serial generation.

A useful additional repository test is an **exhaustive miniature demonstration** over a small field, for example `GF(2^4)` with `k = 4`. There are only `2^16 = 65536` coefficient seeds, so one can enumerate the entire seed space and visibly verify that, at any four selected distinct field points, every possible four-tuple occurs exactly once.

That toy test is not needed for the proof, but it is an unusually transparent executable illustration of the theorem used by the real generator.

---

# 13. Portability and build

The implementation is C++17 and is intended to build with MSVC, GCC, and
Clang.

On x86/x86-64 it detects `PCLMULQDQ` at runtime. GCC/Clang use a
function-specific target attribute plus `__builtin_cpu_supports("pclmul")`,
so the complete program does **not** need to be compiled globally with
`-mpclmul`.

## Direct GCC build

The project is header-only apart from the example/test `main.cpp`, so a direct
Linux build needs only one compilation command:

```bash
g++ -std=c++17 -O3 -pthread main.cpp -o chacha20gf512
./chacha20gf512 1000000 4
```

The second command-line argument is the thread count used by the parallel
benchmark.

To force the portable shift/XOR multiplication path for validation:

```bash
g++ -std=c++17 -O3 -pthread -DCHACHA20GF512_FORCE_PORTABLE \
    main.cpp -o chacha20gf512-portable
```

The portable path computes the same field and the same output stream; it is
only slower.

## CMake build

A small `CMakeLists.txt` is supplied for users who prefer a conventional
cross-platform build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

On multi-configuration generators such as Visual Studio, `--config Release`
selects the optimized configuration. On Linux with Makefiles or Ninja,
`CMAKE_BUILD_TYPE=Release` selects the release flags at configure time.

To build the portable GF multiplication path through CMake:

```bash
cmake -S . -B build-portable \
    -DCMAKE_BUILD_TYPE=Release \
    -DCHACHA20GF512_FORCE_PORTABLE=ON
cmake --build build-portable --config Release
```

CMake links the platform requirements used by the convenience layer:
`Threads::Threads` for the parallel wrapper and `bcrypt` on Windows for
`BCryptGenRandom`.

All seed serialization and `uint64_t` output interpretation are explicitly
little-endian, so reproducible numeric streams do not depend on host
endianness.

---

# 14. API example

```cpp
#include "ChaCha20GF512FFT.h"

uint8_t chacha_seed[32] = { /* 32 independent seed bytes */ };
uint8_t gf_seed[4096]   = { /* full GF coefficient seed */ };

ChaCha20GF512FFT rng(chacha_seed, gf_seed);

uint64_t a = rng.next_int();
uint64_t b = rng.next_int();

rng.seek(1000000);
uint64_t x = rng.next_int();
```

For deterministic tests only:

```cpp
ChaCha20GF512FFT rng(123456789ULL);
```

For practical OS-backed seeding and parallel bulk generation, include the optional convenience headers `ChaCha20GF512Seed.h` and `ChaCha20GF512Parallel.h`. They do not change the mathematical generator.

The short constructor is a convenience API and must not be confused with the full theoretical initialization model.

---

# 15. What is *not* claimed

ChaCha20GF512 does **not** claim:

- to be empirically “more random than ChaCha20”;
- to be the fastest or smallest PRNG;
- full independence of an unlimited output stream;
- exact independence beyond 512 distinct output words;
- exact independence for arbitrary 32768 individual bit positions spread over more than 512 words;
- that the GF component increases ChaCha20's cryptographic key strength;
- that a short seed expanded to 4128 bytes gains 33024 bits of entropy;
- backtracking resistance after compromise of a stateful DRBG state;
- a maximal-period theorem for the combined generator;
- cryptographic audit or suitability as a drop-in system CSPRNG;
- novelty of ChaCha20, polynomial `k`-wise independence, additive FFTs, or XOR-combiner reasoning;
- universal mathematical optimality across every conceivable RNG property.

The phrase **“Maximizing What Can Be Proven About a PRNG”** describes the design motivation, not a theorem that no other construction can optimize some different criterion further.

---

# 16. FAQ / likely objections

## “Isn't this pointless? ChaCha20 already passes statistical tests.”

For most practical applications, ChaCha20 alone is already more than sufficient.

That is exactly why this project asks a different question: **after empirical testing stops providing a useful ranking, can useful exact guarantees still be added?**

ChaCha20GF512 adds an information-theoretic property that a 256-bit-key ChaCha20 family with fixed/public stream ID provably cannot have.

## “Isn't the GF part just linear and therefore weak?”

Yes, deliberately.

The GF polynomial is not intended to hide its structure. Its value is that its bounded-independence property is exact and easy to prove. Beyond 512 words it has algebraic relations, and enough unmasked evaluations determine the polynomial.

ChaCha20 is the component intended to mask that structure computationally.

## “Why not use the GF polynomial alone?”

Because exact bounded independence and computational pseudorandomness are different properties.

GF512 gives an exact theorem through 512 output words, but is algebraically predictable beyond its seed dimension. ChaCha20 provides the complementary computational property.

## “Why not derive the 4096-byte GF seed from the ChaCha key?”

Because a deterministic expansion of a 256-bit key reaches at most `2^256` GF states.

The exact 512-wise theorem requires access to the full `2^32768` coefficient space. Expanding a short seed can produce excellent pseudorandom-looking bytes, but it cannot create the missing information-theoretic seed space.

## “Why no ChaCha rekeying?”

Rekeying is a mechanism used by some stateful DRBGs for properties such as backtracking resistance after state compromise.

ChaCha20GF512 instead uses standard counter-indexed ChaCha20 as a pseudorandom stream component. Backtracking resistance is orthogonal to the exact distribution question studied here and would complicate random access and the construction's analysis without strengthening the stated 512-wise theorem.

## “Can XOR with GF512 weaken ChaCha20?”

Under the independent-seed model used for the computational reduction, an efficient distinguisher for `ChaCha20 XOR GF512` would give an efficient distinguisher for ChaCha20 itself.

The 64-bit convenience constructor does not satisfy that independent-seed model and therefore carries no such formal claim.

## “What if ChaCha20 is broken in the future?”

The computational pseudorandomness claim would need to be reconsidered.

The exact 512-wise independence theorem supplied by the independently seeded GF component would remain true. That is a limited but genuine form of graceful degradation.

## “Why exactly 512?”

`k = 512` was chosen as a practical/theoretical balance:

- 512 independent 64-bit coefficients require exactly 4096 seed bytes;
- 4096 bytes are still small enough to fit comfortably in modern private caches;
- the guarantee covers 32768 bits across any 512 selected output words;
- a 512-point additive FFT maps naturally to a 9-stage power-of-two transform.

Other values of `k` are mathematically possible. The current implementation fixes 512 to keep the code and the claim concrete.

## “Is the idea mathematically new?”

No claim of novelty is made for the ingredients.

Random-polynomial constructions of `k`-wise independent variables are classical. ChaCha20 is established. Additive FFTs over characteristic-two finite fields are established. XOR-combiner arguments are established.

The project is an experimental composition of these ideas around a specific design goal: **combine computational pseudorandomness with a high-order exact distribution guarantee and make both roles explicit in a small implementation.**

---

# 17. Repository layout

A minimal repository can remain deliberately small:

```text
README.md
CMakeLists.txt
ChaCha20GF512FFT.h
ChaCha20GF512Parallel.h
ChaCha20GF512Seed.h
main.cpp
```

The roles are intentionally separated:

- `CMakeLists.txt` provides an optional cross-platform build;
- `ChaCha20GF512FFT.h` contains the actual generator and its mathematical core;
- `ChaCha20GF512Parallel.h` adds optional bulk parallel generation of the same
  logical stream;
- `ChaCha20GF512Seed.h` adds practical Windows/Linux OS-backed seed acquisition;
- `main.cpp` is both example program and validation/benchmark driver.

The core generator has no dependency on threading, OS random APIs, or
`RDRAND`.

## Example program flow

The supplied `main.cpp` demonstrates the project in roughly this order:

```text
1. ChaCha20 reference-vector test
2. GF fast-vs-portable arithmetic test
3. additive-FFT vs. Horner equivalence test
4. random-access / seek test
5. serial-vs-parallel identity test
6. operating-system seed demonstration
7. optional RDRAND demonstration
8. ChaCha20 single-thread benchmark
9. ChaCha20GF512 single-thread benchmark
10. ChaCha20GF512 parallel benchmark
```

`RDRAND` appears only in the example program. It is not part of
ChaCha20GF512's construction, full-seed theorem, or recommended entropy model.
The example exists simply because modern x86 hardware randomness is useful
and interesting to demonstrate separately.

---

# 18. Future work

The core construction is intentionally kept simple. Useful extensions include:

- an exhaustive small-field toy proof (`GF(2^4)`, `k=4`);
- PractRand and TestU01/BigCrush as implementation sanity checks rather than proofs;
- benchmarks on newer x86 and ARM CPUs;
- SIMD/VPCLMUL acceleration of the additive FFT;
- parameterizing the independence order `k`;
- a separately studied **robust computational combiner** variant, e.g. independently seeded ChaCha20 XOR AES-CTR XOR GF512, where computational pseudorandomness could survive failure of either one of the two cryptographic components;
- more theoretical exploration of complementary unconditional test classes, such as small-bias or space-bounded constructions, without complicating the core generator prematurely.

These are intentionally future directions rather than claims of the present implementation.

---

# 19. References

### ChaCha20

- Daniel J. Bernstein, **“ChaCha, a variant of Salsa20”**, 2008.  
  https://cr.yp.to/chacha/chacha-20080120.pdf

- Y. Nir and A. Langley, **RFC 8439: ChaCha20 and Poly1305 for IETF Protocols**.  
  RFC 8439 documents the IETF layout and explicitly notes the original 64-bit-nonce / 64-bit-counter ChaCha layout.  
  https://www.rfc-editor.org/rfc/rfc8439.html

### k-wise independence

- James Aspnes, **“k-wise Independence”**, notes on the standard random-polynomial construction over finite fields.  
  https://www.cs.yale.edu/homes/aspnes/pinewiki/KwiseIndependence.html

### Additive FFT

- Shuhong Gao and Todd Mateer, **“Additive Fast Fourier Transforms over Finite Fields”**, IEEE Transactions on Information Theory 56(12), 2010.  
  https://www.math.clemson.edu/~sgao/papers/GM10.pdf

- Sian-Jheng Lin, Wei-Ho Chung, and Yunghsiang S. Han, **“Novel Polynomial Basis and Its Application to Reed-Solomon Erasure Codes”**.  
  https://arxiv.org/abs/1404.3458

### GF(2^64) polynomial

- N. Madden, **Generalised SIV Internet-Draft**, table of primitive polynomials and reduction constants; lists  
  `x^64 + x^4 + x^3 + x + 1` and `0x1B` for `GF(2^64)`.  
  https://datatracker.ietf.org/doc/html/draft-madden-generalised-siv-00

### GCC runtime dispatch

- GCC documentation, **x86 Built-in Functions** (`__builtin_cpu_supports`).  
  https://gcc.gnu.org/onlinedocs/gcc/x86-Built-in-Functions.html

---

# Status

The current implementation is a first research-quality prototype with independent arithmetic and FFT-vs-Horner checks. Its main purpose is to make the construction, its exact guarantees, its assumptions, and its limitations explicit enough to be independently discussed and challenged.

Feedback, counterexamples, theoretical corrections, implementation review, and comparisons with related constructions are welcome.
