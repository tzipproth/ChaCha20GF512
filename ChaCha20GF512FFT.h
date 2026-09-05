#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>

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
// 128-bit position / GF(2^128) value types
// =============================================================================

struct ChaCha20GF512Position128
{
    uint64_t lo = 0;
    uint64_t hi = 0;

    constexpr ChaCha20GF512Position128() = default;
    constexpr ChaCha20GF512Position128(uint64_t low, uint64_t high = 0)
        : lo(low), hi(high) {}

    constexpr bool operator==(const ChaCha20GF512Position128& o) const
    {
        return lo == o.lo && hi == o.hi;
    }
    constexpr bool operator!=(const ChaCha20GF512Position128& o) const
    {
        return !(*this == o);
    }
    constexpr bool operator<(const ChaCha20GF512Position128& o) const
    {
        return hi < o.hi || (hi == o.hi && lo < o.lo);
    }

    constexpr bool is_max() const
    {
        return lo == ~uint64_t{0} && hi == ~uint64_t{0};
    }

    // Adds a 64-bit amount. Returns true on overflow past 2^128-1.
    bool add_u64(uint64_t v)
    {
        const uint64_t old = lo;
        lo += v;
        if (lo < old) {
            ++hi;
            if (hi == 0)
                return true;
        }
        return false;
    }

    bool increment()
    {
        return add_u64(1);
    }

    constexpr ChaCha20GF512Position128 shr1() const
    {
        return ChaCha20GF512Position128((lo >> 1) | (hi << 63), hi >> 1);
    }

    constexpr ChaCha20GF512Position128 shr3() const
    {
        return ChaCha20GF512Position128((lo >> 3) | (hi << 61), hi >> 3);
    }

    constexpr ChaCha20GF512Position128 aligned_down_512() const
    {
        return ChaCha20GF512Position128(lo & ~uint64_t{511}, hi);
    }

    constexpr std::size_t low9() const
    {
        return static_cast<std::size_t>(lo & uint64_t{511});
    }

    constexpr std::size_t low10() const
    {
        return static_cast<std::size_t>(lo & uint64_t{1023});
    }
};

struct ChaCha20GF512Field128
{
    uint64_t lo = 0;
    uint64_t hi = 0;

    constexpr ChaCha20GF512Field128() = default;
    constexpr ChaCha20GF512Field128(uint64_t low, uint64_t high = 0)
        : lo(low), hi(high) {}

    constexpr bool operator==(const ChaCha20GF512Field128& o) const
    {
        return lo == o.lo && hi == o.hi;
    }
    constexpr bool operator!=(const ChaCha20GF512Field128& o) const
    {
        return !(*this == o);
    }
    constexpr bool is_zero() const
    {
        return lo == 0 && hi == 0;
    }
};

static inline ChaCha20GF512Field128 operator^(ChaCha20GF512Field128 a,
                                                ChaCha20GF512Field128 b)
{
    return ChaCha20GF512Field128(a.lo ^ b.lo, a.hi ^ b.hi);
}
static inline ChaCha20GF512Field128& operator^=(ChaCha20GF512Field128& a,
                                                 ChaCha20GF512Field128 b)
{
    a.lo ^= b.lo;
    a.hi ^= b.hi;
    return a;
}

// =============================================================================
// ChaCha20Counter128
// =============================================================================
//
// A position-indexed ChaCha20 component.  A 64-bit stream_id is first mapped
// to a subkey with HChaCha20 (128-bit HChaCha input = stream_id || 0).  The
// resulting key then uses words 12..15 as one 128-bit block index.  This keeps
// the entire 128-bit position space available while preserving independent
// logical streams through subkey derivation.
//
// HChaCha20 output words are 0,1,2,3,12,13,14,15 after 20 rounds, without the
// ChaCha feed-forward addition.
class ChaCha20Counter128
{
public:
    using Position = ChaCha20GF512Position128;

    static constexpr std::size_t KEY_BYTES = 32;
    static constexpr std::size_t BLOCK_BYTES = 64;
    static constexpr std::size_t WORDS_PER_BLOCK = 8;

    ChaCha20Counter128()
    {
        uint8_t key[KEY_BYTES] = {};
        init(key, 0);
    }

    explicit ChaCha20Counter128(const uint8_t* key32, uint64_t stream_id = 0)
    {
        init(key32, stream_id);
    }

    uint64_t next_int()
    {
        if (buffer_index_ >= WORDS_PER_BLOCK)
            refill();

        return buffer_[buffer_index_++];
    }

    void generate(uint64_t* dst, std::size_t count)
    {
        // Consume an already-buffered partial block first.
        while (count != 0 && buffer_index_ < WORDS_PER_BLOCK) {
            *dst++ = buffer_[buffer_index_++];
            --count;
        }

        // Full blocks can be written directly to the caller's buffer.
        while (count >= WORDS_PER_BLOCK) {
            generate_current_block(dst);
            increment_counter();
            dst += WORDS_PER_BLOCK;
            count -= WORDS_PER_BLOCK;
            buffer_index_ = WORDS_PER_BLOCK;
        }

        if (count != 0) {
            refill();
            while (count-- != 0)
                *dst++ = buffer_[buffer_index_++];
        }
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

    void seek(Position word_position)
    {
        const Position block_counter = word_position.shr3();
        const std::size_t word_in_block = static_cast<std::size_t>(word_position.lo & 7ULL);

        set_counter(block_counter);
        generate_current_block(buffer_);
        increment_counter();
        buffer_index_ = word_in_block;
        current_bits_ = 0;
        bit_index_ = 64;
    }

    void seek(uint64_t word_position)
    {
        seek(Position(word_position));
    }

private:
    uint32_t state_[16]{};
    uint64_t buffer_[WORDS_PER_BLOCK]{};
    std::size_t buffer_index_ = WORDS_PER_BLOCK;
    uint64_t current_bits_ = 0;
    int bit_index_ = 64;

    static uint32_t load32_le(const uint8_t* p)
    {
        return  static_cast<uint32_t>(p[0])
            | (static_cast<uint32_t>(p[1]) << 8)
            | (static_cast<uint32_t>(p[2]) << 16)
            | (static_cast<uint32_t>(p[3]) << 24);
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

    static void rounds20(uint32_t x[16])
    {
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
    }

    static void hchacha20_subkey(const uint8_t* key32, uint64_t stream_id,
                                 uint32_t out_key[8])
    {
        uint32_t x[16] = {
            0x61707865u, 0x3320646eu, 0x79622d32u, 0x6b206574u,
            0,0,0,0,0,0,0,0,
            static_cast<uint32_t>(stream_id),
            static_cast<uint32_t>(stream_id >> 32),
            0, 0
        };
        for (int i = 0; i < 8; ++i)
            x[4 + i] = load32_le(key32 + 4 * i);

        rounds20(x);
        out_key[0] = x[0];
        out_key[1] = x[1];
        out_key[2] = x[2];
        out_key[3] = x[3];
        out_key[4] = x[12];
        out_key[5] = x[13];
        out_key[6] = x[14];
        out_key[7] = x[15];
    }

    void init(const uint8_t* key32, uint64_t stream_id)
    {
        state_[0] = 0x61707865u;
        state_[1] = 0x3320646eu;
        state_[2] = 0x79622d32u;
        state_[3] = 0x6b206574u;

        uint32_t subkey[8];
        hchacha20_subkey(key32, stream_id, subkey);
        for (int i = 0; i < 8; ++i)
            state_[4 + i] = subkey[i];

        state_[12] = state_[13] = state_[14] = state_[15] = 0;
        buffer_index_ = WORDS_PER_BLOCK;
        current_bits_ = 0;
        bit_index_ = 64;
    }

    void set_counter(Position counter)
    {
        state_[12] = static_cast<uint32_t>(counter.lo);
        state_[13] = static_cast<uint32_t>(counter.lo >> 32);
        state_[14] = static_cast<uint32_t>(counter.hi);
        state_[15] = static_cast<uint32_t>(counter.hi >> 32);
        buffer_index_ = WORDS_PER_BLOCK;
    }

    void increment_counter()
    {
        ++state_[12];
        if (state_[12] == 0) {
            ++state_[13];
            if (state_[13] == 0) {
                ++state_[14];
                if (state_[14] == 0)
                    ++state_[15];
            }
        }
    }

    void generate_current_block(uint64_t out[WORDS_PER_BLOCK]) const
    {
        uint32_t x[16];
        std::memcpy(x, state_, sizeof(x));
        rounds20(x);
        for (int i = 0; i < 16; ++i)
            x[i] += state_[i];
        for (std::size_t i = 0; i < WORDS_PER_BLOCK; ++i)
            out[i] = static_cast<uint64_t>(x[2 * i])
                   | (static_cast<uint64_t>(x[2 * i + 1]) << 32);
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
//   Let j = floor(i/2).  Then
//   R(2j)   = ChaCha20(subkey(stream_id), 2j)   XOR low64(P(j))
//   R(2j+1) = ChaCha20(subkey(stream_id), 2j+1) XOR high64(P(j))
//
// P is a degree-511 polynomial over GF(2^128), defined by the irreducible
// polynomial x^128 + x^7 + x^2 + x + 1.  The 512 coefficients require an
// 8192-byte perfectly uniform GF seed.  Any 512 distinct field positions give
// 512 independent uniform GF(2^128) values.  Exposing the low/high 64-bit
// halves as adjacent output words preserves exact 512-wise independence of
// arbitrary uint64_t output positions.
//
// The 512 coefficients are converted once from monomial basis to a normalized
// subspace/novel basis.  Each aligned block of 512 consecutive GF evaluation
// points is evaluated by the 9-stage additive FFT and supplies 1024 uint64_t
// output words by using both 64-bit halves of every field value.
class ChaCha20GF512FFT
{
public:
    using Position = ChaCha20GF512Position128;
    using Field = ChaCha20GF512Field128;

    static constexpr std::size_t CHACHA_SEED_BYTES = 32;
    static constexpr std::size_t GF_COEFFICIENTS = 512;
    static constexpr std::size_t GF_ELEMENT_BYTES = 16;
    static constexpr std::size_t GF_SEED_BYTES = GF_COEFFICIENTS * GF_ELEMENT_BYTES;
    static constexpr std::size_t FULL_SEED_BYTES = CHACHA_SEED_BYTES + GF_SEED_BYTES;
    static constexpr std::size_t FFT_LOGN = 9;
    static constexpr std::size_t FFT_SIZE = std::size_t{1} << FFT_LOGN;
    static constexpr std::size_t OUTPUT_WORDS_PER_FFT = FFT_SIZE * 2;

    ChaCha20GF512FFT(const uint8_t* chacha_seed32,
                    const uint8_t* gf_seed8192,
                    uint64_t stream_id = 0)
        : chacha_(chacha_seed32, stream_id)
    {
        init_gf_from_bytes(gf_seed8192);
        finish_initialization();
    }

    explicit ChaCha20GF512FFT(const uint8_t* full_seed8224,
                             uint64_t stream_id = 0)
        : chacha_(full_seed8224, stream_id)
    {
        init_gf_from_bytes(full_seed8224 + CHACHA_SEED_BYTES);
        finish_initialization();
    }

    // Convenience constructor only. It creates only 2^64 possible complete
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
        chacha_ = ChaCha20Counter128(chacha_seed, stream_id);
        init_gf_from_bytes(gf_seed);
        finish_initialization();
    }

    ChaCha20GF512FFT()
        : ChaCha20GF512FFT(uint64_t{0})
    {
    }

    uint64_t next_int()
    {
        if (exhausted_)
            throw std::overflow_error("ChaCha20GF512FFT: 128-bit word-position space exhausted");

        const Position eval_position = word_position_.shr1();
        ensure_gf_eval_block(eval_position);
        const Field g = gf_block_[eval_position.low9()];
        const uint64_t c = chacha_.next_int();
        const uint64_t out = c ^ ((word_position_.lo & 1ULL) ? g.hi : g.lo);

        if (word_position_.is_max())
            exhausted_ = true;
        else
            word_position_.increment();

        return out;
    }

    void generate(uint64_t* dst, std::size_t count)
    {
        if (count == 0)
            return;
        if (exhausted_)
            throw std::overflow_error("ChaCha20GF512FFT: 128-bit word-position space exhausted");

        // Validate the last requested word before writing anything.
        Position last = word_position_;
        if (last.add_u64(static_cast<uint64_t>(count - 1)))
            throw std::overflow_error("ChaCha20GF512FFT: requested range exceeds 128-bit position space");

        while (count != 0) {
            const Position eval_position = word_position_.shr1();
            ensure_gf_eval_block(eval_position);
            std::size_t eval_offset = eval_position.low9();
            const std::size_t first_half = static_cast<std::size_t>(word_position_.lo & 1ULL);
            const std::size_t available = (FFT_SIZE - eval_offset) * 2 - first_half;
            const std::size_t chunk = std::min<std::size_t>(count, available);

            chacha_.generate(dst, chunk);

            std::size_t i = 0;
            if (first_half != 0 && i < chunk) {
                dst[i++] ^= gf_block_[eval_offset++].hi;
            }
            while (i + 1 < chunk) {
                const Field g = gf_block_[eval_offset++];
                dst[i++] ^= g.lo;
                dst[i++] ^= g.hi;
            }
            if (i < chunk)
                dst[i] ^= gf_block_[eval_offset].lo;

            Position chunk_last = word_position_;
            const bool overflow = chunk_last.add_u64(static_cast<uint64_t>(chunk - 1));
            (void)overflow; // already ruled out by the full-range validation above
            if (chunk_last.is_max()) {
                exhausted_ = true;
                word_position_ = chunk_last;
            }
            else {
                chunk_last.increment();
                word_position_ = chunk_last;
            }

            dst += chunk;
            count -= chunk;
        }
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

    void seek(Position word_position)
    {
        chacha_.seek(word_position);
        word_position_ = word_position;
        exhausted_ = false;
        // Keep a cached FFT block when seeking within the same aligned block.
        // ensure_gf_eval_block() is keyed by gf_block_base_ and will replace it if needed.
        current_bits_ = 0;
        bit_index_ = 64;
    }

    void seek(uint64_t word_position)
    {
        seek(Position(word_position));
    }

    Position position128() const
    {
        return word_position_;
    }

    static constexpr std::size_t gf_seed_state_bytes() { return GF_SEED_BYTES; }
    static constexpr std::size_t gf_cached_block_bytes() { return FFT_SIZE * sizeof(Field); }

    // Slow direct reference evaluation of the full GF(2^128) value. Test use only.
    Field debug_gf_horner128(Position position) const
    {
        return horner_reference(position_to_field(position));
    }

    // FFT result for an arbitrary position, full GF(2^128) value. Test use only.
    Field debug_gf_fft128(Position position)
    {
        ensure_gf_eval_block(position);
        return gf_block_[position.low9()];
    }

    // Backward-friendly low-64 projection helpers for ordinary 64-bit positions.
    uint64_t debug_gf_horner(uint64_t position) const
    {
        return debug_gf_horner128(Position(position)).lo;
    }
    uint64_t debug_gf_fft(uint64_t position)
    {
        return debug_gf_fft128(Position(position)).lo;
    }

    static Field debug_gf_mul_portable(Field a, Field b)
    {
        return gf_mul_portable(a, b);
    }

    static Field debug_gf_mul_fast(Field a, Field b)
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

private:
    // x^128 + x^7 + x^2 + x + 1 -> x^128 == x^7+x^2+x+1.
    static constexpr uint64_t GF_REDUCTION = 0x87ULL;

    ChaCha20Counter128 chacha_;

    // Retained unconditionally so the header-only class has one stable layout
    // in every translation unit.  Each element is now 128 bits.
    Field gf_monomial_[GF_COEFFICIENTS]{};
    Field gf_novel_[GF_COEFFICIENTS]{};
    Field gf_block_[FFT_SIZE]{};
    Position gf_block_base_{};
    bool gf_block_valid_ = false;

    Field subspace_[FFT_LOGN][FFT_LOGN]{};
    Field gamma_[FFT_LOGN]{};
    Field gamma_inv_[FFT_LOGN]{};
    Field psi_beta_[FFT_LOGN][FFT_LOGN]{};

    Position word_position_{};
    uint64_t current_bits_ = 0;
    int bit_index_ = 64;
    bool use_pclmul_ = false;
    bool exhausted_ = false;

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

    void init_gf_from_bytes(const uint8_t* seed8192)
    {
        for (std::size_t i = 0; i < GF_COEFFICIENTS; ++i) {
            const Field v(load64_le(seed8192 + 16 * i),
                          load64_le(seed8192 + 16 * i + 8));
            gf_novel_[i] = v;
            gf_monomial_[i] = v;
        }
    }

    static Field position_to_field(Position p)
    {
        return Field(p.lo, p.hi);
    }

    static Field gf_mul_portable(Field a, Field b)
    {
        Field r;
        for (int i = 0; i < 128; ++i) {
            if (b.lo & 1ULL)
                r ^= a;

            b.lo = (b.lo >> 1) | (b.hi << 63);
            b.hi >>= 1;

            const uint64_t carry = a.hi >> 63;
            a.hi = (a.hi << 1) | (a.lo >> 63);
            a.lo <<= 1;
            if (carry)
                a.lo ^= GF_REDUCTION;
        }
        return r;
    }

#if defined(CHACHA20GF512FFT_MSVC_X86)
    static Field gf_mul_pclmul(Field a, Field b)
    {
        const __m128i va0 = _mm_cvtsi64_si128(static_cast<long long>(a.lo));
        const __m128i va1 = _mm_cvtsi64_si128(static_cast<long long>(a.hi));
        const __m128i vb0 = _mm_cvtsi64_si128(static_cast<long long>(b.lo));
        const __m128i vb1 = _mm_cvtsi64_si128(static_cast<long long>(b.hi));
        const __m128i p00 = _mm_clmulepi64_si128(va0, vb0, 0x00);
        const __m128i p11 = _mm_clmulepi64_si128(va1, vb1, 0x00);
        const __m128i va01 = _mm_cvtsi64_si128(static_cast<long long>(a.lo ^ a.hi));
        const __m128i vb01 = _mm_cvtsi64_si128(static_cast<long long>(b.lo ^ b.hi));
        __m128i px = _mm_clmulepi64_si128(va01, vb01, 0x00);
        px = _mm_xor_si128(px, _mm_xor_si128(p00, p11));
        return reduce_cl_product(p00, px, p11);
    }
#elif defined(CHACHA20GF512FFT_GNU_X86)
    CHACHA20GF512FFT_TARGET_PCLMUL
    static Field gf_mul_pclmul(Field a, Field b)
    {
        const __m128i va0 = _mm_cvtsi64_si128(static_cast<long long>(a.lo));
        const __m128i va1 = _mm_cvtsi64_si128(static_cast<long long>(a.hi));
        const __m128i vb0 = _mm_cvtsi64_si128(static_cast<long long>(b.lo));
        const __m128i vb1 = _mm_cvtsi64_si128(static_cast<long long>(b.hi));
        const __m128i p00 = _mm_clmulepi64_si128(va0, vb0, 0x00);
        const __m128i p11 = _mm_clmulepi64_si128(va1, vb1, 0x00);
        const __m128i va01 = _mm_cvtsi64_si128(static_cast<long long>(a.lo ^ a.hi));
        const __m128i vb01 = _mm_cvtsi64_si128(static_cast<long long>(b.lo ^ b.hi));
        __m128i px = _mm_clmulepi64_si128(va01, vb01, 0x00);
        px = _mm_xor_si128(px, _mm_xor_si128(p00, p11));
        return reduce_cl_product(p00, px, p11);
    }
#endif

#if defined(CHACHA20GF512FFT_MSVC_X86) || defined(CHACHA20GF512FFT_GNU_X86)
#  if defined(CHACHA20GF512FFT_GNU_X86)
    CHACHA20GF512FFT_TARGET_PCLMUL
#  endif
    static Field reduce_cl_product(__m128i p00, __m128i cross, __m128i p11)
    {
        const uint64_t p00lo = static_cast<uint64_t>(_mm_cvtsi128_si64(p00));
        const uint64_t p00hi = static_cast<uint64_t>(_mm_cvtsi128_si64(_mm_srli_si128(p00, 8)));
        const uint64_t xlo   = static_cast<uint64_t>(_mm_cvtsi128_si64(cross));
        const uint64_t xhi   = static_cast<uint64_t>(_mm_cvtsi128_si64(_mm_srli_si128(cross, 8)));
        const uint64_t p11lo = static_cast<uint64_t>(_mm_cvtsi128_si64(p11));
        const uint64_t p11hi = static_cast<uint64_t>(_mm_cvtsi128_si64(_mm_srli_si128(p11, 8)));

        uint64_t c0 = p00lo;
        uint64_t c1 = p00hi ^ xlo;
        const uint64_t c2 = p11lo ^ xhi;
        const uint64_t c3 = p11hi;

        // Fold H*x^128 with x^128 = 1+x+x^2+x^7.  H=(c2,c3).
        uint64_t t0 = 0, t1 = 0, t2 = 0;
        t0 ^= c2; t1 ^= c3; // shift 0

        t0 ^= c2 << 1;
        t1 ^= (c3 << 1) | (c2 >> 63);
        t2 ^= c3 >> 63;

        t0 ^= c2 << 2;
        t1 ^= (c3 << 2) | (c2 >> 62);
        t2 ^= c3 >> 62;

        t0 ^= c2 << 7;
        t1 ^= (c3 << 7) | (c2 >> 57);
        t2 ^= c3 >> 57;

        c0 ^= t0;
        c1 ^= t1;

        // t2 has at most 7 significant bits and represents terms x^(128+j).
        // Fold those once more; the result is below degree 14.
        const uint64_t f = t2 ^ (t2 << 1) ^ (t2 << 2) ^ (t2 << 7);
        c0 ^= f;

        return Field(c0, c1);
    }
#endif

    static bool cpu_has_pclmul()
    {
#if defined(CHACHA20GF512_FORCE_PORTABLE)
        return false;
#else
        // CPUID capability is process-invariant; do not query it in hot test loops.
        static const bool supported = []() {
#if defined(CHACHA20GF512FFT_MSVC_X86)
            int regs[4] = {};
            __cpuid(regs, 1);
            return (static_cast<uint32_t>(regs[2]) & (1u << 1)) != 0;
#elif defined(CHACHA20GF512FFT_GNU_X86)
            __builtin_cpu_init();
            return __builtin_cpu_supports("pclmul") != 0;
#else
            return false;
#endif
        }();
        return supported;
#endif
    }

    Field gf_mul(Field a, Field b) const
    {
#if defined(CHACHA20GF512FFT_MSVC_X86) || defined(CHACHA20GF512FFT_GNU_X86)
        if (use_pclmul_)
            return gf_mul_pclmul(a, b);
#endif
        return gf_mul_portable(a, b);
    }

    // a^(2^128-2).  Initialization-only; simplicity matters more than speed.
    static Field gf_inv_portable(Field a)
    {
        Field r(1, 0);
        Field p = a;
        for (unsigned bit = 0; bit < 128; ++bit) {
            if (bit != 0) // exponent 2^128-2 has bits 1..127 set, bit 0 clear
                r = gf_mul_portable(r, p);
            p = gf_mul_portable(p, p);
        }
        return r;
    }

    Field eval_subspace(std::size_t level, Field x) const
    {
        Field r;
        Field xp = x;
        for (std::size_t j = 0; j <= level; ++j) {
            r ^= gf_mul(subspace_[level][j], xp);
            if (j != level)
                xp = gf_mul(xp, xp);
        }
        return r;
    }

    void build_subspace_polynomials()
    {
        subspace_[0][0] = Field(1, 0); // s_0(X)=X

        for (std::size_t i = 0; i < FFT_LOGN; ++i) {
            const Field beta(uint64_t{1} << i, 0);
            gamma_[i] = eval_subspace(i, beta);
            gamma_inv_[i] = gf_inv_portable(gamma_[i]);

            for (std::size_t j = i + 1; j < FFT_LOGN; ++j) {
                const Field sj = eval_subspace(i, Field(uint64_t{1} << j, 0));
                psi_beta_[i][j] = gf_mul(sj, gamma_inv_[i]);
            }

            if (i + 1 == FFT_LOGN)
                break;

            Field next[FFT_LOGN] = {};
            next[0] = gf_mul(gamma_[i], subspace_[i][0]);
            for (std::size_t j = 1; j <= i; ++j) {
                next[j] = gf_mul(subspace_[i][j - 1], subspace_[i][j - 1])
                        ^ gf_mul(gamma_[i], subspace_[i][j]);
            }
            next[i + 1] = gf_mul(subspace_[i][i], subspace_[i][i]);
            for (std::size_t j = 0; j <= i + 1; ++j)
                subspace_[i + 1][j] = next[j];
        }
    }

    void monomial_to_novel(Field* data, std::size_t logn)
    {
        if (logn == 0)
            return;

        const std::size_t level = logn - 1;
        const std::size_t half = std::size_t{1} << level;
        const std::size_t n = half * 2;

        for (std::size_t idx = n; idx-- > half;) {
            const Field q = data[idx];
            const std::size_t qi = idx - half;
            for (std::size_t j = 0; j < level; ++j)
                data[qi + (std::size_t{1} << j)] ^=
                    gf_mul(q, subspace_[level][j]);
            data[idx] = gf_mul(gamma_[level], q);
        }

        monomial_to_novel(data, logn - 1);
        monomial_to_novel(data + half, logn - 1);
    }

    void finish_initialization()
    {
        use_pclmul_ = cpu_has_pclmul();
        build_subspace_polynomials();
        monomial_to_novel(gf_novel_, FFT_LOGN);
        word_position_ = Position();
        gf_block_base_ = Position();
        gf_block_valid_ = false;
        current_bits_ = 0;
        bit_index_ = 64;
        exhausted_ = false;
    }

    Field horner_reference(Field x) const
    {
        Field y = gf_monomial_[GF_COEFFICIENTS - 1];
        for (std::size_t i = GF_COEFFICIENTS - 1; i-- > 0;)
            y = gf_mul_portable(y, x) ^ gf_monomial_[i];
        return y;
    }

    void ensure_gf_eval_block(Position position)
    {
        const Position base = position.aligned_down_512();
        if (gf_block_valid_ && gf_block_base_ == base)
            return;

        std::memcpy(gf_block_, gf_novel_, sizeof(gf_block_));
        fft_block(base);
        gf_block_base_ = base;
        gf_block_valid_ = true;
    }

    void prepare_base_t(Position base, Field base_t[FFT_LOGN]) const
    {
        const Field x = position_to_field(base);
        for (std::size_t i = 0; i < FFT_LOGN; ++i) {
            const Field s = eval_subspace(i, x);
            base_t[i] = gf_mul(s, gamma_inv_[i]);
        }
    }

    void fft_recursive_portable(Field* data,
                                std::size_t logn,
                                uint16_t high_mask,
                                const Field base_t[FFT_LOGN]) const
    {
        if (logn == 0)
            return;

        const std::size_t level = logn - 1;
        const std::size_t half = std::size_t{1} << level;

        Field t = base_t[level];
        for (std::size_t j = level + 1; j < FFT_LOGN; ++j) {
            if (high_mask & (uint16_t(1) << j))
                t ^= psi_beta_[level][j];
        }

        for (std::size_t i = 0; i < half; ++i) {
            const Field upper = data[i + half];
            const Field left = data[i] ^ gf_mul_portable(t, upper);
            data[i] = left;
            data[i + half] = left ^ upper;
        }

        fft_recursive_portable(data, logn - 1, high_mask, base_t);
        fft_recursive_portable(data + half, logn - 1,
                               static_cast<uint16_t>(high_mask | (uint16_t(1) << level)),
                               base_t);
    }

#if defined(CHACHA20GF512FFT_MSVC_X86) || defined(CHACHA20GF512FFT_GNU_X86)
#  if defined(CHACHA20GF512FFT_GNU_X86)
    CHACHA20GF512FFT_TARGET_PCLMUL
#  endif
    void fft_recursive_pclmul(Field* data,
                              std::size_t logn,
                              uint16_t high_mask,
                              const Field base_t[FFT_LOGN]) const
    {
        if (logn == 0)
            return;

        const std::size_t level = logn - 1;
        const std::size_t half = std::size_t{1} << level;

        Field t = base_t[level];
        for (std::size_t j = level + 1; j < FFT_LOGN; ++j) {
            if (high_mask & (uint16_t(1) << j))
                t ^= psi_beta_[level][j];
        }

        for (std::size_t i = 0; i < half; ++i) {
            const Field upper = data[i + half];
            const Field left = data[i] ^ gf_mul_pclmul(t, upper);
            data[i] = left;
            data[i + half] = left ^ upper;
        }

        fft_recursive_pclmul(data, logn - 1, high_mask, base_t);
        fft_recursive_pclmul(data + half, logn - 1,
                             static_cast<uint16_t>(high_mask | (uint16_t(1) << level)),
                             base_t);
    }
#endif

    void fft_block(Position base)
    {
        Field base_t[FFT_LOGN];
        prepare_base_t(base, base_t);
#if defined(CHACHA20GF512FFT_MSVC_X86) || defined(CHACHA20GF512FFT_GNU_X86)
        if (use_pclmul_) {
            fft_recursive_pclmul(gf_block_, FFT_LOGN, 0, base_t);
            return;
        }
#endif
        fft_recursive_portable(gf_block_, FFT_LOGN, 0, base_t);
    }
};

#undef CHACHA20GF512FFT_TARGET_PCLMUL
#undef CHACHA20GF512FFT_GNU_X86
#undef CHACHA20GF512FFT_MSVC_X86
