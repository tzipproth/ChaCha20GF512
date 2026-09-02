#include "ChaCha20GF512FFT.h"
#include "ChaCha20GF512Parallel.h"
#include "ChaCha20GF512Seed.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

static bool test_chacha_zero_vector()
{
    static const uint64_t expected[8] = {
        0x903df1a0ade0b876ULL, 0x28bd8653e56a5d40ULL,
        0x1aed8da0b819d2bdULL, 0xc70d778bccef36a8ULL,
        0x8d4857517c5941daULL, 0x374ad8b83fe02477ULL,
        0x1ca11815f4b8436aULL, 0x8665eeb269b687c3ULL
    };
    uint8_t key[32] = {};
    ChaCha20Counter64 c(key, 0);
    for (int i = 0; i < 8; ++i)
        if (c.next_int() != expected[i]) return false;
    return true;
}

static bool test_gf_arithmetic()
{
    uint64_t s = 0xA0761D6478BD642FULL;
    for (int i = 0; i < 200000; ++i) {
        const uint64_t a = sm64(s);
        const uint64_t b = sm64(s);
        const uint64_t c = sm64(s);
        const uint64_t p = ChaCha20GF512FFT::debug_gf_mul_portable(a, b);
        const uint64_t f = ChaCha20GF512FFT::debug_gf_mul_fast(a, b);
        if (p != f) return false;
        if (ChaCha20GF512FFT::debug_gf_mul_fast(a, b ^ c) !=
            (ChaCha20GF512FFT::debug_gf_mul_fast(a, b) ^
             ChaCha20GF512FFT::debug_gf_mul_fast(a, c))) return false;
        if (ChaCha20GF512FFT::debug_gf_mul_fast(a, 1) != a) return false;
    }
    return true;
}

static bool test_fft_vs_horner()
{
    const uint64_t fixed[] = {
        0,1,2,7,8,255,256,510,511,512,513,1023,1024,
        1234,65535,65536,0x123456789ABC0000ULL,
        0xFFFFFFFFFFFFFE00ULL,0xFFFFFFFFFFFFFFFFULL
    };

    for (int seed_no = 0; seed_no < 4; ++seed_no) {
        uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES];
        uint64_t ss = 0x123456789ABCDEF0ULL ^
                      (0x9E3779B97F4A7C15ULL * uint64_t(seed_no + 1));
        for (std::size_t i = 0; i < sizeof(seed); i += 8)
            store64_le(seed + i, sm64(ss));
        ChaCha20GF512FFT g(seed);

        for (uint64_t p : fixed)
            if (g.debug_gf_fft(p) != g.debug_gf_horner(p)) return false;

        uint64_t r = 0xD1B54A32D192ED03ULL ^ uint64_t(seed_no);
        for (int i = 0; i < 500; ++i) {
            const uint64_t p = sm64(r);
            if (g.debug_gf_fft(p) != g.debug_gf_horner(p)) return false;
        }
    }
    return true;
}

static bool test_seek()
{
    uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES];
    make_reproducible_full_seed(seed);
    ChaCha20GF512FFT seq(seed), rnd(seed);
    constexpr uint64_t POS = 12345;
    uint64_t expected = 0;
    for (uint64_t i = 0; i <= POS; ++i) expected = seq.next_int();
    rnd.seek(POS);
    return expected == rnd.next_int();
}

static bool test_parallel_identity(unsigned threads)
{
    uint8_t seed[ChaCha20GF512FFT::FULL_SEED_BYTES];
    make_reproducible_full_seed(seed);

    constexpr uint64_t START = 123;
    constexpr std::size_t COUNT = 20000 + 37;
    std::vector<uint64_t> serial(COUNT), parallel(COUNT);

    ChaCha20GF512FFT s(seed);
    s.seek(START);
    for (std::size_t i = 0; i < COUNT; ++i)
        serial[i] = s.next_int();

    ChaCha20GF512Parallel p(seed);
    p.fill_parallel_at(parallel.data(), parallel.size(), START, threads);
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
// NOT claim 33024 bits of fresh physical entropy. It is kept out of the RNG and
// seed helper APIs so the core design is not coupled to x86 hardware behavior.
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
    for (int i = 0; i < 4; ++i) store64_le(key + i*8, sm64(s));
    ChaCha20Counter64 rng(key);
    volatile uint64_t sink = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (uint64_t i = 0; i < count; ++i) sink ^= rng.next_int();
    const auto t1 = std::chrono::steady_clock::now();
    const double sec = std::chrono::duration<double>(t1-t0).count();
    std::printf("%-32s %9.3f M uint64/s  %9.3f MiB/s  sink=%016llx\n",
        "ChaCha20Counter64", count/sec/1e6, count*8.0/sec/(1024*1024),
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
    const double sec = std::chrono::duration<double>(t1-t0).count();
    std::printf("%-32s %9.3f M uint64/s  %9.3f MiB/s  sink=%016llx\n",
        "ChaCha20GF512 FFT single", count/sec/1e6, count*8.0/sec/(1024*1024),
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
    const double sec = std::chrono::duration<double>(t1-t0).count();
    char name[64];
    std::snprintf(name, sizeof(name), "ChaCha20GF512 FFT parallel x%u", threads);
    std::printf("%-32s %9.3f M uint64/s  %9.3f MiB/s  sink=%016llx\n",
        name, count/sec/1e6, count*8.0/sec/(1024*1024),
        (unsigned long long)sink);
}

int main(int argc, char** argv)
{
    uint64_t count = 1000000ULL;
    unsigned threads = 4;
    if (argc > 1) count = std::strtoull(argv[1], nullptr, 10);
    if (argc > 2) threads = static_cast<unsigned>(std::strtoul(argv[2], nullptr, 10));

    std::printf("PCLMUL runtime support:    %s\n",
                ChaCha20GF512FFT::debug_cpu_has_pclmul() ? "yes" : "no");
    std::printf("sizeof(ChaCha20GF512FFT):    %zu bytes\n",
                sizeof(ChaCha20GF512FFT));
    std::printf("ChaCha 64/64 test vector:  %s\n",
                test_chacha_zero_vector() ? "OK" : "FAILED");
    std::printf("GF arithmetic fast/ref:    %s\n",
                test_gf_arithmetic() ? "OK" : "FAILED");
    std::printf("FFT vs Horner:             %s\n",
                test_fft_vs_horner() ? "OK" : "FAILED");
    std::printf("Hybrid seek test:          %s\n",
                test_seek() ? "OK" : "FAILED");
    std::printf("Parallel identity x%u:      %s\n", threads,
                test_parallel_identity(threads) ? "OK" : "FAILED");

    demo_seeding();

    std::printf("Iterations: %llu\n", (unsigned long long)count);
    bench_chacha(count);
    bench_single(count);
    bench_parallel(count, threads);
    return 0;
}
