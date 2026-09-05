#include "ChaCha20GF512FFT.h"
#include "ChaCha20GF512Parallel.h"
#include "ChaCha20GF512Seed.h"
#include "toy_kwise_test.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
#  include <immintrin.h>
#  include <intrin.h>
#  define DEMO_MSVC_X86 1
#elif (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#  include <immintrin.h>
#  define DEMO_GNU_X86 1
#  define DEMO_TARGET_RDRND __attribute__((target("rdrnd")))
#else
#  define DEMO_TARGET_RDRND
#endif

using Position128 = ChaCha20GF512FFT::Position;
using Field128 = ChaCha20GF512FFT::Field;

static uint64_t sm64(uint64_t& s)
{
    uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static void store64_le(uint8_t* p, uint64_t v)
{
    for (int i = 0; i < 8; ++i)
        p[i] = static_cast<uint8_t>(v >> (8 * i));
}

static void make_reproducible_full_seed(
    uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES])
{
    uint64_t s = 0x123456789ABCDEF0ULL;
    for (std::size_t i = 0; i < ChaCha20GF512FFT::FULL_SEED_BYTES; i += 8)
        store64_le(seed + i, sm64(s));
}

static Position128 add_pos(Position128 p, uint64_t delta)
{
    if (p.add_u64(delta))
        throw std::overflow_error("test position overflow");
    return p;
}

// Regression vector for the new ChaCha20Counter128 layout:
// zero 256-bit key, stream_id = 0, block position = 0.
// The stream-id subkey is HChaCha20(key, stream_id || 0), then words 12..15
// are used as one 128-bit block index.
static bool test_chacha_zero_vector()
{
    static const uint64_t expected[8] = {
        0xd1013fbf182ad0bcULL, 0xacfda8a730de9292ULL,
        0xc72c00a6505eb6a4ULL, 0xd5c31ac9f7d2d62cULL,
        0xcfbfd2aae0838f72ULL, 0xddae8fb52d2dbd9aULL,
        0x139bc03fd85d0165ULL, 0x0f8e9e014310271eULL
    };

    uint8_t key[32] = {};
    ChaCha20Counter128 c(key, 0);
    for (int i = 0; i < 8; ++i)
        if (c.next_int() != expected[i]) return false;
    return true;
}

static bool test_hybrid_regression_vector()
{
    static const uint64_t expected[16] = {
        0xde9aff20072a8ee4ULL, 0x680b3c09a5cb1a1bULL,
        0xcf4adae3ecb51785ULL, 0xa0e824c25c835030ULL,
        0x12a3a8e92708f081ULL, 0x94b5d72863c6bc56ULL,
        0xf7362ba3c5e9b8c7ULL, 0xaa1bdb38b38755e3ULL,
        0xb148bb6c24e2eeddULL, 0x98a32053d2550c82ULL,
        0x1086db2d3e2075a1ULL, 0xe94c0e238000c700ULL,
        0xe0faed830bd798a1ULL, 0x43303620284f904cULL,
        0xfbd4a81b79200de2ULL, 0x04b4c74ded93e16bULL
    };

    uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES];
    make_reproducible_full_seed(seed);
    ChaCha20GF512FFT rng(seed);
    for (uint64_t v : expected)
        if (rng.next_int() != v) return false;
    return true;
}

static bool test_gf_arithmetic()
{
    uint64_t s = 0xA0761D6478BD642FULL;
    const Field128 one{1, 0};

    for (int i = 0; i < 100000; ++i) {
        const Field128 a{sm64(s), sm64(s)};
        const Field128 b{sm64(s), sm64(s)};
        const Field128 c{sm64(s), sm64(s)};

        const Field128 p = ChaCha20GF512FFT::debug_gf_mul_portable(a, b);
        const Field128 f = ChaCha20GF512FFT::debug_gf_mul_fast(a, b);
        if (p != f) return false;

        if (ChaCha20GF512FFT::debug_gf_mul_fast(a, b ^ c) !=
            (ChaCha20GF512FFT::debug_gf_mul_fast(a, b) ^
             ChaCha20GF512FFT::debug_gf_mul_fast(a, c))) return false;

        if (ChaCha20GF512FFT::debug_gf_mul_fast(a, one) != a) return false;
    }
    return true;
}

static bool test_fft_vs_horner()
{
    const Position128 fixed[] = {
        {0,0}, {1,0}, {2,0}, {7,0}, {8,0}, {255,0}, {256,0},
        {510,0}, {511,0}, {512,0}, {513,0}, {1023,0}, {1024,0},
        {1234,0}, {65535,0}, {65536,0},
        {0x123456789ABC0000ULL,0},
        {0xFFFFFFFFFFFFFE00ULL,0},
        {0xFFFFFFFFFFFFFFFFULL,0},
        {0,1},
        {1,1},
        {0x123456789ABCDEF0ULL,0x0123456789ABCDEFULL},
        {0xFFFFFFFFFFFFFE00ULL,0xFEDCBA9876543210ULL}
    };

    for (int seed_no = 0; seed_no < 4; ++seed_no) {
        uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES];
        uint64_t ss = 0x123456789ABCDEF0ULL ^
                      (0x9E3779B97F4A7C15ULL * uint64_t(seed_no + 1));
        for (std::size_t i = 0; i < sizeof(seed); i += 8)
            store64_le(seed + i, sm64(ss));
        ChaCha20GF512FFT g(seed);

        for (const Position128 p : fixed)
            if (g.debug_gf_fft128(p) != g.debug_gf_horner128(p)) return false;

        uint64_t r = 0xD1B54A32D192ED03ULL ^ uint64_t(seed_no);
        for (int i = 0; i < 500; ++i) {
            const Position128 p{sm64(r), sm64(r)};
            if (g.debug_gf_fft128(p) != g.debug_gf_horner128(p)) return false;
        }
    }
    return true;
}

static bool test_seek()
{
    uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES];
    make_reproducible_full_seed(seed);

    // Ordinary backward-compatible 64-bit seek check.
    {
        ChaCha20GF512FFT seq(seed), rnd(seed);
        constexpr uint64_t POS = 12345;
        uint64_t expected = 0;
        for (uint64_t i = 0; i <= POS; ++i) expected = seq.next_int();
        rnd.seek(POS);
        if (expected != rnd.next_int()) return false;
    }

    // New 128-bit check: walk across the low-64-bit wrap and compare every
    // sequential output with an independent direct seek to the same position.
    {
        const Position128 start{std::numeric_limits<uint64_t>::max() - 12ULL,
                                0x0123456789ABCDEFULL};
        ChaCha20GF512FFT seq(seed), direct(seed);
        seq.seek(start);
        for (uint64_t i = 0; i < 40; ++i) {
            const uint64_t a = seq.next_int();
            const Position128 p = add_pos(start, i);
            direct.seek(p);
            const uint64_t b = direct.next_int();
            if (a != b) return false;
        }
    }

    return true;
}

static bool test_bulk_identity()
{
    uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES];
    make_reproducible_full_seed(seed);

    struct Case { Position128 start; std::size_t count; };
    const Case cases[] = {
        {{0, 0}, 4099},
        {{1, 0}, 2051},
        {{1023, 0}, 3077},
        {{1024, 0}, 2049},
        {{std::numeric_limits<uint64_t>::max() - 1500ULL, 0x123456789ABCDEF0ULL}, 4096},
        {{std::numeric_limits<uint64_t>::max() - 100ULL,
          std::numeric_limits<uint64_t>::max()}, 101}
    };

    for (const Case& tc : cases) {
        std::vector<uint64_t> scalar(tc.count), bulk(tc.count);
        ChaCha20GF512FFT a(seed), b(seed);
        a.seek(tc.start);
        b.seek(tc.start);
        for (std::size_t i = 0; i < tc.count; ++i)
            scalar[i] = a.next_int();
        b.generate(bulk.data(), bulk.size());
        if (scalar != bulk) return false;
    }

    // A one-word request at the final representable position is valid, then
    // both APIs must report exhaustion rather than wrapping to zero.
    const Position128 last{std::numeric_limits<uint64_t>::max(),
                           std::numeric_limits<uint64_t>::max()};
    ChaCha20GF512FFT a(seed), b(seed);
    a.seek(last);
    b.seek(last);
    uint64_t bv = 0;
    b.generate(&bv, 1);
    if (bv != a.next_int()) return false;

    bool a_threw = false, b_threw = false;
    try { (void)a.next_int(); } catch (const std::overflow_error&) { a_threw = true; }
    try { b.generate(&bv, 1); } catch (const std::overflow_error&) { b_threw = true; }
    return a_threw && b_threw;
}

static bool test_parallel_identity(unsigned threads)
{
    uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES];
    make_reproducible_full_seed(seed);

    constexpr std::size_t COUNT = 20000 + 37;
    std::vector<uint64_t> serial(COUNT), parallel(COUNT);

    // Deliberately cross the 64-bit low-half boundary. This verifies that the
    // parallel wrapper preserves the same logical 128-bit-position stream.
    const Position128 start{std::numeric_limits<uint64_t>::max() - 700ULL,
                            0x1111222233334444ULL};

    ChaCha20GF512FFT s(seed);
    s.seek(start);
    for (std::size_t i = 0; i < COUNT; ++i)
        serial[i] = s.next_int();

    ChaCha20GF512Parallel p(seed);
    p.fill_parallel_at(parallel.data(), parallel.size(), start, threads);
    return serial == parallel;
}

static bool cpu_has_rdrand()
{
#if defined(DEMO_MSVC_X86)
    int regs[4] = {};
    __cpuid(regs, 1);
    return (static_cast<uint32_t>(regs[2]) & (1u << 30)) != 0;
#elif defined(DEMO_GNU_X86)
    __builtin_cpu_init();
    return __builtin_cpu_supports("rdrnd") != 0;
#else
    return false;
#endif
}

#if defined(DEMO_MSVC_X86)
static bool rdrand64(uint64_t& out)
{
#  if defined(_M_X64)
    return _rdrand64_step(reinterpret_cast<unsigned __int64*>(&out)) != 0;
#  else
    unsigned int lo = 0, hi = 0;
    if (!_rdrand32_step(&lo) || !_rdrand32_step(&hi)) return false;
    out = uint64_t(lo) | (uint64_t(hi) << 32);
    return true;
#  endif
}
#elif defined(DEMO_GNU_X86)
DEMO_TARGET_RDRND
static bool rdrand64(uint64_t& out)
{
#  if defined(__x86_64__)
    unsigned long long v = 0;
    const int ok = _rdrand64_step(&v);
    out = static_cast<uint64_t>(v);
    return ok != 0;
#  else
    unsigned int lo = 0, hi = 0;
    if (!_rdrand32_step(&lo) || !_rdrand32_step(&hi)) return false;
    out = uint64_t(lo) | (uint64_t(hi) << 32);
    return true;
#  endif
}
#else
static bool rdrand64(uint64_t&) { return false; }
#endif

// Demo only. RDRAND is the output of a hardware DRBG; this deliberately does
// NOT claim FULL_SEED_BYTES*8 bits of fresh physical entropy. It is kept out
// of the RNG and seed helper APIs so the core design is not coupled to x86.
static bool make_rdrand_demo_seed(
    uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES])
{
    if (!cpu_has_rdrand()) return false;
    for (std::size_t i = 0; i < ChaCha20GF512FFT::FULL_SEED_BYTES; i += 8) {
        uint64_t v = 0;
        bool ok = false;
        for (int retry = 0; retry < 16 && !ok; ++retry)
            ok = rdrand64(v);
        if (!ok) return false;
        store64_le(seed + i, v);
    }
    return true;
}

static void demo_seeding()
{
    try {
        const auto os_seed = chacha20gf512_seed::make_os_full_seed();
        ChaCha20GF512FFT os_rng(os_seed.data());
        std::printf("OS random seed demo:       OK  first=%016llx\n",
                    (unsigned long long)os_rng.next_int());
    }
    catch (...) {
        std::printf("OS random seed demo:       FAILED\n");
    }

    uint8_t hw_seed[ChaCha20GF512FFT::FULL_SEED_BYTES];
    if (make_rdrand_demo_seed(hw_seed)) {
        ChaCha20GF512FFT hw_rng(hw_seed);
        std::printf("RDRAND seed demo:          OK  first=%016llx\n",
                    (unsigned long long)hw_rng.next_int());
    }
    else {
        std::printf("RDRAND seed demo:          unavailable/failed\n");
    }
}

static void bench_chacha(uint64_t count)
{
    uint8_t key[32] = {};
    uint64_t s = 0x123456789ABCDEF0ULL;
    for (int i = 0; i < 4; ++i) store64_le(key + i * 8, sm64(s));
    ChaCha20Counter128 rng(key);
    volatile uint64_t sink = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (uint64_t i = 0; i < count; ++i) sink ^= rng.next_int();
    const auto t1 = std::chrono::steady_clock::now();
    const double sec = std::chrono::duration<double>(t1 - t0).count();
    const double dcount = static_cast<double>(count);
    std::printf("%-32s %9.3f M uint64/s  %9.3f MiB/s  sink=%016llx\n",
        "ChaCha20Counter128", dcount / sec / 1e6, dcount * 8.0 / sec / (1024 * 1024),
        (unsigned long long)sink);
}

static void bench_single(uint64_t count)
{
    uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES];
    make_reproducible_full_seed(seed);
    ChaCha20GF512FFT rng(seed);
    volatile uint64_t sink = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (uint64_t i = 0; i < count; ++i) sink ^= rng.next_int();
    const auto t1 = std::chrono::steady_clock::now();
    const double sec = std::chrono::duration<double>(t1 - t0).count();
    const double dcount = static_cast<double>(count);
    std::printf("%-32s %9.3f M uint64/s  %9.3f MiB/s  sink=%016llx\n",
        "ChaCha20GF512 FFT single", dcount / sec / 1e6, dcount * 8.0 / sec / (1024 * 1024),
        (unsigned long long)sink);
}

static void bench_bulk(uint64_t count)
{
    uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES];
    make_reproducible_full_seed(seed);
    ChaCha20GF512FFT rng(seed);
    std::vector<uint64_t> out(static_cast<std::size_t>(count));

    const auto t0 = std::chrono::steady_clock::now();
    rng.generate(out.data(), out.size());
    const auto t1 = std::chrono::steady_clock::now();

    uint64_t sink = 0;
    for (uint64_t v : out) sink ^= v;
    const double sec = std::chrono::duration<double>(t1 - t0).count();
    const double dcount = static_cast<double>(count);
    std::printf("%-32s %9.3f M uint64/s  %9.3f MiB/s  sink=%016llx\n",
        "ChaCha20GF512 FFT bulk", dcount / sec / 1e6, dcount * 8.0 / sec / (1024 * 1024),
        (unsigned long long)sink);
}

static void bench_parallel(uint64_t count, unsigned threads)
{
    uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES];
    make_reproducible_full_seed(seed);
    ChaCha20GF512Parallel rng(seed);
    std::vector<uint64_t> out(static_cast<std::size_t>(count));

    const auto t0 = std::chrono::steady_clock::now();
    rng.fill_parallel(out.data(), out.size(), threads);
    const auto t1 = std::chrono::steady_clock::now();

    uint64_t sink = 0;
    for (uint64_t v : out) sink ^= v; // outside timed region
    const double sec = std::chrono::duration<double>(t1 - t0).count();
    char name[64];
    std::snprintf(name, sizeof(name), "ChaCha20GF512 FFT parallel x%u", threads);
    const double dcount = static_cast<double>(count);
    std::printf("%-32s %9.3f M uint64/s  %9.3f MiB/s  sink=%016llx\n",
        name, dcount / sec / 1e6, dcount * 8.0 / sec / (1024 * 1024),
        (unsigned long long)sink);
}

int main(int argc, char** argv)
{
    uint64_t count = 1000000ULL;
    unsigned threads = 4;
    if (argc > 1) count = std::strtoull(argv[1], nullptr, 10);
    if (argc > 2) threads = static_cast<unsigned>(std::strtoul(argv[2], nullptr, 10));
    const bool full_toy = argc > 3 && std::strcmp(argv[3], "full") == 0;

    std::printf("PCLMUL runtime support:    %s\n",
        ChaCha20GF512FFT::debug_cpu_has_pclmul() ? "yes" : "no");
    std::printf("sizeof(ChaCha20GF512FFT):  %zu bytes\n",
        sizeof(ChaCha20GF512FFT));
    std::printf("Full seed size:            %zu bytes\n",
        ChaCha20GF512FFT::FULL_SEED_BYTES);
    std::printf("GF seed size:              %zu bytes\n",
        ChaCha20GF512FFT::GF_SEED_BYTES);
    std::printf("ChaCha128 regression vec:  %s\n",
        test_chacha_zero_vector() ? "OK" : "FAILED");
    std::printf("Hybrid regression vector: %s\n",
        test_hybrid_regression_vector() ? "OK" : "FAILED");
    std::printf("GF(2^128) fast/ref:        %s\n",
        test_gf_arithmetic() ? "OK" : "FAILED");
    std::printf("GF128 FFT vs Horner:       %s\n",
        test_fft_vs_horner() ? "OK" : "FAILED");
    std::printf("Hybrid 128-bit seek:       %s\n",
        test_seek() ? "OK" : "FAILED");
    std::printf("Bulk vs next_int:          %s\n",
        test_bulk_identity() ? "OK" : "FAILED");
    std::printf("Parallel identity x%u:      %s\n", threads,
        test_parallel_identity(threads) ? "OK" : "FAILED");

    demo_seeding();

    std::printf("Iterations: %llu\n", (unsigned long long)count);
    bench_chacha(count);
    bench_single(count);
    bench_bulk(count);
    bench_parallel(count, threads);

    // ------------------------------------------------------------------
    // Exhaustive toy-model verification of the exact k-wise claim.
    // The production construction now uses GF(2^128) with a 64-bit output
    // projection. The toy models exhaust smaller fields and verify:
    // (A) exact k-wise independence,
    // (B) sharpness at k+1,
    // (C) invariance under a fixed XOR mask,
    // (D) a reducible-polynomial negative control,
    // (E) that projecting field values to fewer output bits preserves k-wise
    //     independence with exactly the expected multiplicities.
    // ------------------------------------------------------------------
    std::printf("\nToy-model k-wise verification (%s mode, exhaustive seed-space enumeration,\n"
        "GF(2^4) k=4%s):\n",
        full_toy ? "full" : "quick",
        full_toy ? " and GF(2^8) k=3" : "");
    std::fflush(stdout);
    const bool toy_ok = run_toy_kwise_test(full_toy, /*gf8_sets=*/4, stdout);
    std::printf("Toy-model k-wise check:    %s\n", toy_ok ? "OK" : "FAILED");

    return toy_ok ? 0 : 1;
}
