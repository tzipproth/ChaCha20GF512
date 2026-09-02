#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>

#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_AMD64) || defined(_M_IX86))
#  include <intrin.h>
#  include <wmmintrin.h>
#  define CHACHA20GF512FFT_MSVC_X86 1
#elif (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#  include <wmmintrin.h>
#  define CHACHA20GF512FFT_GNU_X86 1
#  define CHACHA20GF512FFT_TARGET_PCLMUL __attribute__((target("pclmul,sse2")))
#else
#  define CHACHA20GF512FFT_TARGET_PCLMUL
#endif

// =============================================================================
// ChaCha20Counter64
// =============================================================================
// Original ChaCha layout:
//   256-bit key + 64-bit block counter + 64-bit stream id / nonce.
//
// This is deliberately a counter stream, not a state-updating DRBG.  It gives
// a simple position-indexed ChaCha20 component with random access.
class ChaCha20Counter64
{
public:
    static constexpr std::size_t KEY_BYTES = 32;
    static constexpr std::size_t BLOCK_BYTES = 64;
    static constexpr std::size_t WORDS_PER_BLOCK = 8;

    ChaCha20Counter64()
    {
        uint8_t key[KEY_BYTES] = {};
        init(key, 0);
    }

    explicit ChaCha20Counter64(const uint8_t* key32, uint64_t stream_id = 0)
    {
        init(key32, stream_id);
    }

    uint64_t next_int()
    {
        if (buffer_index_ > BLOCK_BYTES - sizeof(uint64_t))
            refill();

        const uint64_t x = load64_le(buffer_ + buffer_index_);
        buffer_index_ += sizeof(uint64_t);
        return x;
    }

    int next_bit()
    {
        if (bit_index_ == 64) {
            current_bits_ = next_int();
            bit_index_ = 0;
        }
        const int bit = static_cast<int>((current_bits_ >> bit_index_) & 1ULL);
        ++bit_index_;
        return bit;
    }

    void seek(uint64_t word_position)
    {
        const uint64_t block_counter = word_position / WORDS_PER_BLOCK;
        const std::size_t word_in_block = static_cast<std::size_t>(
            word_position % WORDS_PER_BLOCK);

        set_counter(block_counter);
        generate_current_block(buffer_);
        increment_counter();
        buffer_index_ = word_in_block * sizeof(uint64_t);
        current_bits_ = 0;
        bit_index_ = 64;
    }

private:
    uint32_t state_[16]{};
    uint8_t buffer_[BLOCK_BYTES]{};
    std::size_t buffer_index_ = BLOCK_BYTES;
    uint64_t current_bits_ = 0;
    int bit_index_ = 64;

    static uint32_t load32_le(const uint8_t* p)
    {
        return  static_cast<uint32_t>(p[0])
            | (static_cast<uint32_t>(p[1]) << 8)
            | (static_cast<uint32_t>(p[2]) << 16)
            | (static_cast<uint32_t>(p[3]) << 24);
    }

    static uint64_t load64_le(const uint8_t* p)
    {
        return  static_cast<uint64_t>(p[0])
            | (static_cast<uint64_t>(p[1]) << 8)
            | (static_cast<uint64_t>(p[2]) << 16)
            | (static_cast<uint64_t>(p[3]) << 24)
            | (static_cast<uint64_t>(p[4]) << 32)
            | (static_cast<uint64_t>(p[5]) << 40)
            | (static_cast<uint64_t>(p[6]) << 48)
            | (static_cast<uint64_t>(p[7]) << 56);
    }

    static void store32_le(uint8_t* p, uint32_t v)
    {
        p[0] = static_cast<uint8_t>(v);
        p[1] = static_cast<uint8_t>(v >> 8);
        p[2] = static_cast<uint8_t>(v >> 16);
        p[3] = static_cast<uint8_t>(v >> 24);
    }

    static uint32_t rotl32(uint32_t x, int n)
    {
        return (x << n) | (x >> (32 - n));
    }

    static void quarterround(uint32_t& a, uint32_t& b,
                             uint32_t& c, uint32_t& d)
    {
        a += b; d ^= a; d = rotl32(d, 16);
        c += d; b ^= c; b = rotl32(b, 12);
        a += b; d ^= a; d = rotl32(d, 8);
        c += d; b ^= c; b = rotl32(b, 7);
    }

    void init(const uint8_t* key32, uint64_t stream_id)
    {
        state_[0] = 0x61707865;
        state_[1] = 0x3320646e;
        state_[2] = 0x79622d32;
        state_[3] = 0x6b206574;

        for (int i = 0; i < 8; ++i)
            state_[4 + i] = load32_le(key32 + 4 * i);

        state_[12] = 0;
        state_[13] = 0;
        state_[14] = static_cast<uint32_t>(stream_id);
        state_[15] = static_cast<uint32_t>(stream_id >> 32);

        buffer_index_ = BLOCK_BYTES;
        current_bits_ = 0;
        bit_index_ = 64;
    }

    void set_counter(uint64_t counter)
    {
        state_[12] = static_cast<uint32_t>(counter);
        state_[13] = static_cast<uint32_t>(counter >> 32);
        buffer_index_ = BLOCK_BYTES;
    }

    void increment_counter()
    {
        ++state_[12];
        if (state_[12] == 0)
            ++state_[13];
    }

    void generate_current_block(uint8_t out[BLOCK_BYTES]) const
    {
        uint32_t x[16];
        std::memcpy(x, state_, sizeof(x));

        for (int i = 0; i < 10; ++i) {
            quarterround(x[0], x[4], x[8],  x[12]);
            quarterround(x[1], x[5], x[9],  x[13]);
            quarterround(x[2], x[6], x[10], x[14]);
            quarterround(x[3], x[7], x[11], x[15]);

            quarterround(x[0], x[5], x[10], x[15]);
            quarterround(x[1], x[6], x[11], x[12]);
            quarterround(x[2], x[7], x[8],  x[13]);
            quarterround(x[3], x[4], x[9],  x[14]);
        }

        for (int i = 0; i < 16; ++i)
            x[i] += state_[i];
        for (int i = 0; i < 16; ++i)
            store32_le(out + 4 * i, x[i]);
    }

    void refill()
    {
        generate_current_block(buffer_);
        increment_counter();
        buffer_index_ = 0;
    }
};

// =============================================================================
// ChaCha20GF512FFT
// =============================================================================
//
//   R(i) = ChaCha20(key, stream_id, i) XOR P(i)
//
// P is a degree-511 polynomial over GF(2^64).  The 512 monomial coefficients
// supplied by the 4096-byte GF seed are converted ONCE during initialization
// into a normalized subspace / novel basis.  Thereafter each aligned block of
// 512 consecutive field points is evaluated by a 9-stage additive FFT.
//
// The point set of an aligned block B = 512*k is
//
//   { B, B+1, ..., B+511 }
// = { B XOR j : 0 <= j < 512 }
// = B + span_F2(1,2,4,...,256),
//
// an affine 9-dimensional F2-subspace of GF(2^64).  This is exactly the domain
// on which the additive FFT operates.
//
// The mathematical generator is UNCHANGED relative to direct Horner
// evaluation.  The FFT is only a faster way of computing the same P(i).
//
// Runtime multiplication count per 512-word GF block:
//   9 FFT stages * 256 butterflies = 2304 GF multiplications,
//   plus 90 GF multiplications to evaluate the nine affine base constants.
// Total: 2394 / 512 ~= 4.68 GF multiplications per output word, versus 511
// with direct Horner evaluation.  The monomial->novel conversion is paid once
// during initialization and does not affect steady-state throughput.
class ChaCha20GF512FFT
{
public:
    static constexpr std::size_t CHACHA_SEED_BYTES = 32;
    static constexpr std::size_t GF_COEFFICIENTS = 512;
    static constexpr std::size_t GF_SEED_BYTES = GF_COEFFICIENTS * sizeof(uint64_t);
    static constexpr std::size_t FULL_SEED_BYTES = CHACHA_SEED_BYTES + GF_SEED_BYTES;
    static constexpr std::size_t FFT_LOGN = 9;
    static constexpr std::size_t FFT_SIZE = 1u << FFT_LOGN;

    ChaCha20GF512FFT(const uint8_t* chacha_seed32,
                   const uint8_t* gf_seed4096,
                   uint64_t stream_id = 0)
        : chacha_(chacha_seed32, stream_id)
    {
        init_gf_from_bytes(gf_seed4096);
        finish_initialization();
    }

    explicit ChaCha20GF512FFT(const uint8_t* full_seed4128,
                            uint64_t stream_id = 0)
        : chacha_(full_seed4128, stream_id)
    {
        init_gf_from_bytes(full_seed4128 + CHACHA_SEED_BYTES);
        finish_initialization();
    }

    // Convenience constructor only.  It creates only 2^64 possible complete
    // initial states and therefore DOES NOT provide the full-seed-space exact
    // 512-wise guarantee or 256 bits of ChaCha key entropy.
    explicit ChaCha20GF512FFT(uint64_t seed, uint64_t stream_id = 0)
        : chacha_()
    {
        uint8_t chacha_seed[CHACHA_SEED_BYTES];
        uint8_t gf_seed[GF_SEED_BYTES];
        uint64_t sm = seed;
        fill_from_splitmix64(chacha_seed, sizeof(chacha_seed), sm);
        fill_from_splitmix64(gf_seed, sizeof(gf_seed), sm);
        chacha_ = ChaCha20Counter64(chacha_seed, stream_id);
        init_gf_from_bytes(gf_seed);
        finish_initialization();
    }

    ChaCha20GF512FFT()
        : ChaCha20GF512FFT(uint64_t{0})
    {
    }

    uint64_t next_int()
    {
        ensure_gf_block(word_position_);
        const uint64_t g = gf_block_[static_cast<std::size_t>(word_position_ & (FFT_SIZE - 1))];
        const uint64_t c = chacha_.next_int();
        ++word_position_;
        return c ^ g;
    }

    int next_bit()
    {
        if (bit_index_ == 64) {
            current_bits_ = next_int();
            bit_index_ = 0;
        }
        const int bit = static_cast<int>((current_bits_ >> bit_index_) & 1ULL);
        ++bit_index_;
        return bit;
    }

    void seek(uint64_t word_position)
    {
        chacha_.seek(word_position);
        word_position_ = word_position;
        gf_block_valid_ = false;
        current_bits_ = 0;
        bit_index_ = 64;
    }

    static constexpr std::size_t gf_seed_state_bytes() { return GF_SEED_BYTES; }
    static constexpr std::size_t gf_cached_block_bytes() { return FFT_SIZE * sizeof(uint64_t); }

#ifdef CHACHA20GF512_ENABLE_TEST_API
    // Slow direct reference evaluation of P(position).  Test use only.
    uint64_t debug_gf_horner(uint64_t position) const
    {
        return horner_reference(position);
    }

    // FFT result for an arbitrary position.  This may fill one 512-word cache.
    uint64_t debug_gf_fft(uint64_t position)
    {
        ensure_gf_block(position);
        return gf_block_[static_cast<std::size_t>(position & (FFT_SIZE - 1))];
    }

    static uint64_t debug_gf_mul_portable(uint64_t a, uint64_t b)
    {
        return gf_mul_portable(a, b);
    }

    static uint64_t debug_gf_mul_fast(uint64_t a, uint64_t b)
    {
#if defined(CHACHA20GF512FFT_MSVC_X86) || defined(CHACHA20GF512FFT_GNU_X86)
        if (cpu_has_pclmul())
            return gf_mul_pclmul(a, b);
#endif
        return gf_mul_portable(a, b);
    }

    static bool debug_cpu_has_pclmul()
    {
        return cpu_has_pclmul();
    }
#endif

private:
    static constexpr uint64_t GF_REDUCTION = 0x1BULL;
    static constexpr uint64_t INVALID_BLOCK_BASE = ~uint64_t{0};

    static constexpr uint64_t RED_OVER[16] = {
        0x00, 0x1B, 0x2D, 0x36, 0x5A, 0x41, 0x77, 0x6C,
        0xAF, 0xB4, 0x82, 0x99, 0xF5, 0xEE, 0xD8, 0xC3
    };

    ChaCha20Counter64 chacha_;

#ifdef CHACHA20GF512_ENABLE_TEST_API
    // Retained only in test builds so direct Horner evaluation can verify the
    // FFT bit-for-bit.  Production builds do not pay these extra 4096 bytes.
    alignas(64) uint64_t gf_monomial_[GF_COEFFICIENTS]{};
#endif

    // The 4096-byte GF coefficient state, represented in the normalized
    // subspace/novel basis after a one-time invertible initialization transform.
    alignas(64) uint64_t gf_novel_[GF_COEFFICIENTS]{};

    // Cached evaluations of one aligned affine subspace (512 output words).
    alignas(64) uint64_t gf_block_[FFT_SIZE]{};
    uint64_t gf_block_base_ = INVALID_BLOCK_BASE;
    bool gf_block_valid_ = false;

    // Sparse subspace polynomials:
    // s_i(X) = product_{v in V_i} (X-v)
    //        = sum_{j=0..i} subspace_[i][j] * X^(2^j)
    // where V_i = span(1,2,...,2^(i-1)).
    uint64_t subspace_[FFT_LOGN][FFT_LOGN]{};
    uint64_t gamma_[FFT_LOGN]{};       // gamma_i = s_i(beta_i)
    uint64_t gamma_inv_[FFT_LOGN]{};

    // psi_i(beta_j), where psi_i = s_i / s_i(beta_i).  Only j>i is used
    // by the FFT; psi_i(beta_i)=1 by construction.
    uint64_t psi_beta_[FFT_LOGN][FFT_LOGN]{};

    uint64_t word_position_ = 0;
    uint64_t current_bits_ = 0;
    int bit_index_ = 64;
    bool use_pclmul_ = false;

    static uint64_t load64_le(const uint8_t* p)
    {
        return  static_cast<uint64_t>(p[0])
            | (static_cast<uint64_t>(p[1]) << 8)
            | (static_cast<uint64_t>(p[2]) << 16)
            | (static_cast<uint64_t>(p[3]) << 24)
            | (static_cast<uint64_t>(p[4]) << 32)
            | (static_cast<uint64_t>(p[5]) << 40)
            | (static_cast<uint64_t>(p[6]) << 48)
            | (static_cast<uint64_t>(p[7]) << 56);
    }

    static void store64_le(uint8_t* p, uint64_t v)
    {
        for (int i = 0; i < 8; ++i)
            p[i] = static_cast<uint8_t>(v >> (8 * i));
    }

    static uint64_t splitmix64_step(uint64_t& state)
    {
        uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    static void fill_from_splitmix64(uint8_t* dst, std::size_t bytes,
                                      uint64_t& sm_state)
    {
        while (bytes >= 8) {
            const uint64_t v = splitmix64_step(sm_state);
            store64_le(dst, v);
            dst += 8;
            bytes -= 8;
        }
        if (bytes != 0) {
            uint8_t tmp[8];
            store64_le(tmp, splitmix64_step(sm_state));
            std::memcpy(dst, tmp, bytes);
        }
    }

    void init_gf_from_bytes(const uint8_t* seed4096)
    {
        for (std::size_t i = 0; i < GF_COEFFICIENTS; ++i) {
            const uint64_t v = load64_le(seed4096 + 8 * i);
            gf_novel_[i] = v;
#ifdef CHACHA20GF512_ENABLE_TEST_API
            gf_monomial_[i] = v;
#endif
        }
    }

    static inline uint64_t gf_mul_portable(uint64_t a, uint64_t b)
    {
        uint64_t r = 0;
        for (int i = 0; i < 64; ++i) {
            if (b & 1ULL)
                r ^= a;
            b >>= 1;
            const uint64_t carry = a >> 63;
            a <<= 1;
            if (carry)
                a ^= GF_REDUCTION;
        }
        return r;
    }

#if defined(CHACHA20GF512FFT_MSVC_X86)
    static inline uint64_t gf_mul_pclmul(uint64_t a, uint64_t b)
    {
        const __m128i va = _mm_cvtsi64_si128(static_cast<long long>(a));
        const __m128i vb = _mm_cvtsi64_si128(static_cast<long long>(b));
        const __m128i p = _mm_clmulepi64_si128(va, vb, 0x00);
        const uint64_t lo = static_cast<uint64_t>(_mm_cvtsi128_si64(p));
        const uint64_t hi = static_cast<uint64_t>(
            _mm_cvtsi128_si64(_mm_srli_si128(p, 8)));
        return lo ^ hi ^ (hi << 1) ^ (hi << 3) ^ (hi << 4)
                  ^ RED_OVER[hi >> 60];
    }
#elif defined(CHACHA20GF512FFT_GNU_X86)
    CHACHA20GF512FFT_TARGET_PCLMUL
    static inline uint64_t gf_mul_pclmul(uint64_t a, uint64_t b)
    {
        const __m128i va = _mm_cvtsi64_si128(static_cast<long long>(a));
        const __m128i vb = _mm_cvtsi64_si128(static_cast<long long>(b));
        const __m128i p = _mm_clmulepi64_si128(va, vb, 0x00);
        const uint64_t lo = static_cast<uint64_t>(_mm_cvtsi128_si64(p));
        const uint64_t hi = static_cast<uint64_t>(
            _mm_cvtsi128_si64(_mm_srli_si128(p, 8)));
        return lo ^ hi ^ (hi << 1) ^ (hi << 3) ^ (hi << 4)
                  ^ RED_OVER[hi >> 60];
    }
#endif

    static bool cpu_has_pclmul()
    {
#if defined(CHACHA20GF512_FORCE_PORTABLE)
        return false;
#elif defined(CHACHA20GF512FFT_MSVC_X86)
        int regs[4] = {};
        __cpuid(regs, 1);
        return (static_cast<uint32_t>(regs[2]) & (1u << 1)) != 0;
#elif defined(CHACHA20GF512FFT_GNU_X86)
        // GCC/Clang perform runtime dispatch; no global -mpclmul is required.
        __builtin_cpu_init();
        return __builtin_cpu_supports("pclmul") != 0;
#else
        return false;
#endif
    }

    static uint64_t gf_pow_portable(uint64_t a, uint64_t e)
    {
        uint64_t r = 1;
        while (e != 0) {
            if (e & 1ULL)
                r = gf_mul_portable(r, a);
            a = gf_mul_portable(a, a);
            e >>= 1;
        }
        return r;
    }

    static uint64_t gf_inv_portable(uint64_t a)
    {
        // a^(2^64-2).  gamma_i is never zero because beta_i is outside V_i.
        return gf_pow_portable(a, ~uint64_t{1});
    }

    uint64_t eval_subspace_portable(std::size_t level, uint64_t x) const
    {
        uint64_t r = 0;
        uint64_t xp = x;
        for (std::size_t j = 0; j <= level; ++j) {
            r ^= gf_mul_portable(subspace_[level][j], xp);
            if (j != level)
                xp = gf_mul_portable(xp, xp);
        }
        return r;
    }

#if defined(CHACHA20GF512FFT_MSVC_X86)
    uint64_t eval_subspace_pclmul(std::size_t level, uint64_t x) const
    {
        uint64_t r = 0;
        uint64_t xp = x;
        for (std::size_t j = 0; j <= level; ++j) {
            r ^= gf_mul_pclmul(subspace_[level][j], xp);
            if (j != level)
                xp = gf_mul_pclmul(xp, xp);
        }
        return r;
    }
#elif defined(CHACHA20GF512FFT_GNU_X86)
    CHACHA20GF512FFT_TARGET_PCLMUL
    uint64_t eval_subspace_pclmul(std::size_t level, uint64_t x) const
    {
        uint64_t r = 0;
        uint64_t xp = x;
        for (std::size_t j = 0; j <= level; ++j) {
            r ^= gf_mul_pclmul(subspace_[level][j], xp);
            if (j != level)
                xp = gf_mul_pclmul(xp, xp);
        }
        return r;
    }
#endif

    void build_subspace_polynomials()
    {
        // s_0(X) = X.
        subspace_[0][0] = 1;

        for (std::size_t i = 0; i < FFT_LOGN; ++i) {
            const uint64_t beta = uint64_t{1} << i;
            gamma_[i] = eval_subspace_portable(i, beta);
            gamma_inv_[i] = gf_inv_portable(gamma_[i]);

            // Values needed by the FFT's affine-coset butterflies.
            for (std::size_t j = i + 1; j < FFT_LOGN; ++j) {
                const uint64_t sj = eval_subspace_portable(i, uint64_t{1} << j);
                psi_beta_[i][j] = gf_mul_portable(sj, gamma_inv_[i]);
            }

            if (i + 1 == FFT_LOGN)
                break;

            // s_{i+1}(X) = s_i(X)^2 + gamma_i * s_i(X).
            uint64_t next[FFT_LOGN] = {};
            next[0] = gf_mul_portable(gamma_[i], subspace_[i][0]);
            for (std::size_t j = 1; j <= i; ++j) {
                next[j] = gf_mul_portable(subspace_[i][j - 1], subspace_[i][j - 1])
                        ^ gf_mul_portable(gamma_[i], subspace_[i][j]);
            }
            next[i + 1] = gf_mul_portable(subspace_[i][i], subspace_[i][i]);
            for (std::size_t j = 0; j <= i + 1; ++j)
                subspace_[i + 1][j] = next[j];
        }
    }

    // Convert a monomial-basis polynomial segment of length 2^logn in-place
    // into the normalized novel basis generated by psi_i = s_i/gamma_i.
    // This is an initialization-only operation, so the portable multiplier is
    // intentionally used here for simple, architecture-independent behavior.
    void monomial_to_novel(uint64_t* data, std::size_t logn)
    {
        if (logn == 0)
            return;

        const std::size_t level = logn - 1;
        const std::size_t half = std::size_t{1} << level;
        const std::size_t n = half * 2;

        // Long-divide by the monic s_level.  If q is the quotient and r the
        // remainder, then f = r + s*q = r + psi*(gamma*q), so gamma*q is the
        // upper-half coefficient vector in the normalized novel basis.
        for (std::size_t idx = n; idx-- > half;) {
            const uint64_t q = data[idx];
            const std::size_t qi = idx - half;
            for (std::size_t j = 0; j < level; ++j)
                data[qi + (std::size_t{1} << j)] ^=
                    gf_mul_portable(q, subspace_[level][j]);
            data[idx] = gf_mul_portable(gamma_[level], q);
        }

        monomial_to_novel(data, logn - 1);
        monomial_to_novel(data + half, logn - 1);
    }

    void finish_initialization()
    {
        use_pclmul_ = cpu_has_pclmul();
        build_subspace_polynomials();
        monomial_to_novel(gf_novel_, FFT_LOGN);
        word_position_ = 0;
        gf_block_base_ = INVALID_BLOCK_BASE;
        gf_block_valid_ = false;
        current_bits_ = 0;
        bit_index_ = 64;
    }

#ifdef CHACHA20GF512_ENABLE_TEST_API
    uint64_t horner_reference(uint64_t x) const
    {
        uint64_t y = gf_monomial_[GF_COEFFICIENTS - 1];
        for (std::size_t i = GF_COEFFICIENTS - 1; i-- > 0;)
            y = gf_mul_portable(y, x) ^ gf_monomial_[i];
        return y;
    }
#endif

    void ensure_gf_block(uint64_t position)
    {
        const uint64_t base = position & ~uint64_t(FFT_SIZE - 1);
        if (gf_block_valid_ && gf_block_base_ == base)
            return;

        std::memcpy(gf_block_, gf_novel_, sizeof(gf_block_));

#if defined(CHACHA20GF512FFT_MSVC_X86) || defined(CHACHA20GF512FFT_GNU_X86)
        if (use_pclmul_)
            fft_block_pclmul(base);
        else
            fft_block_portable(base);
#else
        fft_block_portable(base);
#endif

        gf_block_base_ = base;
        gf_block_valid_ = true;
    }

    void prepare_base_t_portable(uint64_t base, uint64_t base_t[FFT_LOGN]) const
    {
        for (std::size_t i = 0; i < FFT_LOGN; ++i) {
            const uint64_t s = eval_subspace_portable(i, base);
            base_t[i] = gf_mul_portable(s, gamma_inv_[i]);
        }
    }

    void fft_recursive_portable(uint64_t* data,
                                std::size_t logn,
                                uint16_t high_mask,
                                const uint64_t base_t[FFT_LOGN]) const
    {
        if (logn == 0)
            return;

        const std::size_t level = logn - 1;
        const std::size_t half = std::size_t{1} << level;

        uint64_t t = base_t[level];
        for (std::size_t j = level + 1; j < FFT_LOGN; ++j) {
            if (high_mask & (uint16_t(1) << j))
                t ^= psi_beta_[level][j];
        }

        for (std::size_t i = 0; i < half; ++i) {
            const uint64_t upper = data[i + half];
            const uint64_t left = data[i] ^ gf_mul_portable(t, upper);
            data[i] = left;
            data[i + half] = left ^ upper; // psi_level(beta_level) == 1
        }

        fft_recursive_portable(data, logn - 1, high_mask, base_t);
        fft_recursive_portable(data + half, logn - 1,
                               static_cast<uint16_t>(high_mask | (uint16_t(1) << level)),
                               base_t);
    }

    void fft_block_portable(uint64_t base)
    {
        uint64_t base_t[FFT_LOGN];
        prepare_base_t_portable(base, base_t);
        fft_recursive_portable(gf_block_, FFT_LOGN, 0, base_t);
    }

#if defined(CHACHA20GF512FFT_MSVC_X86)
    void prepare_base_t_pclmul(uint64_t base, uint64_t base_t[FFT_LOGN]) const
    {
        for (std::size_t i = 0; i < FFT_LOGN; ++i) {
            const uint64_t s = eval_subspace_pclmul(i, base);
            base_t[i] = gf_mul_pclmul(s, gamma_inv_[i]);
        }
    }

    void fft_recursive_pclmul(uint64_t* data,
                              std::size_t logn,
                              uint16_t high_mask,
                              const uint64_t base_t[FFT_LOGN]) const
    {
        if (logn == 0)
            return;
        const std::size_t level = logn - 1;
        const std::size_t half = std::size_t{1} << level;
        uint64_t t = base_t[level];
        for (std::size_t j = level + 1; j < FFT_LOGN; ++j)
            if (high_mask & (uint16_t(1) << j))
                t ^= psi_beta_[level][j];
        for (std::size_t i = 0; i < half; ++i) {
            const uint64_t upper = data[i + half];
            const uint64_t left = data[i] ^ gf_mul_pclmul(t, upper);
            data[i] = left;
            data[i + half] = left ^ upper;
        }
        fft_recursive_pclmul(data, logn - 1, high_mask, base_t);
        fft_recursive_pclmul(data + half, logn - 1,
                             static_cast<uint16_t>(high_mask | (uint16_t(1) << level)),
                             base_t);
    }

    void fft_block_pclmul(uint64_t base)
    {
        uint64_t base_t[FFT_LOGN];
        prepare_base_t_pclmul(base, base_t);
        fft_recursive_pclmul(gf_block_, FFT_LOGN, 0, base_t);
    }
#elif defined(CHACHA20GF512FFT_GNU_X86)
    CHACHA20GF512FFT_TARGET_PCLMUL
    void prepare_base_t_pclmul(uint64_t base, uint64_t base_t[FFT_LOGN]) const
    {
        for (std::size_t i = 0; i < FFT_LOGN; ++i) {
            const uint64_t s = eval_subspace_pclmul(i, base);
            base_t[i] = gf_mul_pclmul(s, gamma_inv_[i]);
        }
    }

    CHACHA20GF512FFT_TARGET_PCLMUL
    void fft_recursive_pclmul(uint64_t* data,
                              std::size_t logn,
                              uint16_t high_mask,
                              const uint64_t base_t[FFT_LOGN]) const
    {
        if (logn == 0)
            return;
        const std::size_t level = logn - 1;
        const std::size_t half = std::size_t{1} << level;
        uint64_t t = base_t[level];
        for (std::size_t j = level + 1; j < FFT_LOGN; ++j)
            if (high_mask & (uint16_t(1) << j))
                t ^= psi_beta_[level][j];
        for (std::size_t i = 0; i < half; ++i) {
            const uint64_t upper = data[i + half];
            const uint64_t left = data[i] ^ gf_mul_pclmul(t, upper);
            data[i] = left;
            data[i + half] = left ^ upper;
        }
        fft_recursive_pclmul(data, logn - 1, high_mask, base_t);
        fft_recursive_pclmul(data + half, logn - 1,
                             static_cast<uint16_t>(high_mask | (uint16_t(1) << level)),
                             base_t);
    }

    CHACHA20GF512FFT_TARGET_PCLMUL
    void fft_block_pclmul(uint64_t base)
    {
        uint64_t base_t[FFT_LOGN];
        prepare_base_t_pclmul(base, base_t);
        fft_recursive_pclmul(gf_block_, FFT_LOGN, 0, base_t);
    }
#endif
};

#undef CHACHA20GF512FFT_TARGET_PCLMUL
#undef CHACHA20GF512FFT_GNU_X86
#undef CHACHA20GF512FFT_MSVC_X86
