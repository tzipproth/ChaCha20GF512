#pragma once
// toy_kwise_test.h
//
// Exhaustive brute-force verification of the exact k-wise independence claim
// behind ChaCha20GF512, scaled down to fields whose ENTIRE seed space can be
// enumerated.  The construction is identical to the production generator:
// output at position x is P(x) with x interpreted as a field element, over a
// field defined by a low-weight irreducible polynomial.
//
//   Production: GF(2^64), 512 coefficients, seed space 2^32768 (not enumerable)
//   Toy:        GF(2^4),  4 coefficients,   seed space 2^16
//               GF(2^8),  3 coefficients,   seed space 2^24
//
// Checks, all by exhaustive counting over the seed space:
//   (A) k-wise:    for k distinct positions every k-tuple occurs EXACTLY once
//   (B) sharpness: for k+1 positions only 2^(mk) of 2^(m(k+1)) tuples are reachable
//   (C) masking:   XOR with an arbitrary fixed position-dependent mask (the
//                  role of ChaCha20 in the hybrid) leaves (A) unchanged
//   (D) control:   with a REDUCIBLE polynomial (x^4+1) check (A) must FAIL,
//                  proving the test can detect a broken construction
//
// Modes:
//   quick  (~1 s)  GF(2^4): (A) over ALL 1820 position sets, (B)/(C) sampled,
//                  (D) negative control.  Suitable as a default self-test.
//   full  (~15 s)  additionally (B)/(C) over ALL sets and GF(2^8) with
//                  `gf8_sets` sampled position sets.
//
// This does not test the 64-bit generator; it tests the theorem and the
// construction.  The proof is field-independent, so a defect would show up
// here.  See README, "Toy-model verification".

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

// ---------------------------------------------------------------------------
// Small binary field GF(2^M) with reduction polynomial x^M + LOW.
// ---------------------------------------------------------------------------
template <unsigned M, unsigned LOW>
struct SmallField
{
    static constexpr unsigned SIZE = 1u << M;
    static constexpr unsigned MASK = SIZE - 1;

    static unsigned mul_slow(unsigned a, unsigned b)
    {
        unsigned r = 0;
        for (unsigned i = 0; i < M; ++i) {
            if (b & 1u) r ^= a;
            b >>= 1;
            const unsigned carry = a & (1u << (M - 1));
            a = (a << 1) & MASK;
            if (carry) a ^= LOW;
        }
        return r;
    }

    // Full multiplication table (256x256 bytes for GF(2^8)).
    static const std::vector<uint8_t>& table()
    {
        static const std::vector<uint8_t> t = [] {
            std::vector<uint8_t> v(std::size_t(SIZE) * SIZE);
            for (unsigned a = 0; a < SIZE; ++a)
                for (unsigned b = 0; b < SIZE; ++b)
                    v[std::size_t(a) * SIZE + b] = uint8_t(mul_slow(a, b));
            return v;
        }();
        return t;
    }

    static inline unsigned mul(unsigned a, unsigned b)
    {
        return table()[std::size_t(a) * SIZE + b];
    }

    // Horner evaluation of P(x) = c[0] + c[1] x + ... + c[k-1] x^(k-1).
    static inline unsigned eval(const unsigned* c, unsigned k, unsigned x)
    {
        unsigned y = c[k - 1];
        for (unsigned i = k - 1; i-- > 0;)
            y = mul(y, x) ^ c[i];
        return y;
    }

    // Sanity: the field must actually be a field (multiplicative group of
    // order 2^M - 1).  Checks a^(2^M-1) == 1 for every nonzero a.
    static bool is_field()
    {
        for (unsigned a = 1; a < SIZE; ++a) {
            unsigned r = 1;
            for (unsigned e = 0; e < SIZE - 1; ++e) r = mul_slow(r, a);
            if (r != 1) return false;
        }
        return true;
    }
};

// Decode seed index -> k coefficients (M bits each).
template <unsigned M>
static inline void seed_to_coeffs(uint64_t seed, unsigned k, unsigned* c)
{
    for (unsigned i = 0; i < k; ++i)
        c[i] = unsigned((seed >> (M * i)) & ((1u << M) - 1));
}

// Arbitrary fixed mask standing in for the ChaCha20 stream: any function of
// the position works for the argument, so a cheap integer hash is enough.
static inline unsigned demo_mask(unsigned position, unsigned m_bits)
{
    uint64_t z = (uint64_t(position) + 1) * 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return unsigned((z ^ (z >> 31)) & ((1u << m_bits) - 1));
}

// ---------------------------------------------------------------------------
// (A)/(C): every k-tuple exactly once over the full seed space.
// ---------------------------------------------------------------------------
template <class F, unsigned M>
static bool check_exactly_once(const unsigned* pos, unsigned k, bool masked)
{
    const uint64_t seeds = 1ull << (M * k);
    std::vector<uint8_t> count(std::size_t(seeds), 0); // tuple index -> count
    unsigned c[16];

    for (uint64_t s = 0; s < seeds; ++s) {
        seed_to_coeffs<M>(s, k, c);
        uint64_t tuple = 0;
        for (unsigned i = 0; i < k; ++i) {
            unsigned v = F::eval(c, k, pos[i]);
            if (masked) v ^= demo_mask(pos[i], M);
            tuple |= uint64_t(v) << (M * i);
        }
        if (++count[std::size_t(tuple)] != 1) return false; // seen twice
    }
    // seeds == number of tuples and every tuple hit at most once -> exactly once
    return true;
}

// ---------------------------------------------------------------------------
// (B): number of distinct reachable (k+1)-tuples.
// ---------------------------------------------------------------------------
template <class F, unsigned M>
static uint64_t count_reachable(const unsigned* pos, unsigned k1)
{
    const unsigned k = k1 - 1;
    const uint64_t seeds = 1ull << (M * k);
    const uint64_t tuples = 1ull << (M * k1);
    std::vector<uint8_t> seen(std::size_t((tuples + 7) / 8), 0);
    uint64_t distinct = 0;
    unsigned c[16];

    for (uint64_t s = 0; s < seeds; ++s) {
        seed_to_coeffs<M>(s, k, c);
        uint64_t tuple = 0;
        for (unsigned i = 0; i < k1; ++i)
            tuple |= uint64_t(F::eval(c, k, pos[i])) << (M * i);
        uint8_t& byte = seen[std::size_t(tuple >> 3)];
        const uint8_t bit = uint8_t(1u << (tuple & 7));
        if (!(byte & bit)) { byte |= bit; ++distinct; }
    }
    return distinct;
}

// Enumerate all k-subsets of {0..n-1}.
static bool next_combination(unsigned* idx, unsigned k, unsigned n)
{
    int i = int(k) - 1;
    while (i >= 0 && idx[i] == n - k + unsigned(i)) --i;
    if (i < 0) return false;
    ++idx[i];
    for (unsigned j = unsigned(i) + 1; j < k; ++j) idx[j] = idx[j - 1] + 1;
    return true;
}

static double seconds_since(std::chrono::steady_clock::time_point t0)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}


namespace toy_kwise_detail {

// GF(2^4), k = 4
template <unsigned LOW>
static bool gf16_check_A(unsigned* pos_out_fail)
{
    using F = SmallField<4, LOW>;
    unsigned pos[4] = {0, 1, 2, 3};
    do {
        if (!check_exactly_once<F, 4>(pos, 4, false)) {
            if (pos_out_fail) std::copy(pos, pos + 4, pos_out_fail);
            return false;
        }
    } while (next_combination(pos, 4, F::SIZE));
    return true;
}

} // namespace toy_kwise_detail

// Returns true if all checks pass.  `out` may be nullptr for silent operation.
static bool run_toy_kwise_test(bool full, unsigned gf8_sets, std::FILE* out)
{
    using namespace toy_kwise_detail;
    auto say = [out](const char* fmt, auto... args) {
        if (!out) return;
        if constexpr (sizeof...(args) == 0) std::fputs(fmt, out);
        else std::fprintf(out, fmt, args...);
    };
    const auto t0 = std::chrono::steady_clock::now();
    bool ok = true;

    {
        using F = SmallField<4, 0x3>;
        constexpr unsigned M = 4, K = 4;
        say("  GF(2^4) k=%u seed space 2^%u%s\n", K, M * K, F::is_field() ? "" : "  FIELD CHECK FAILED");
        ok &= F::is_field();

        unsigned failpos[4];
        if (gf16_check_A<0x3>(failpos))
            say("    (A) k-wise:    all 1820 position sets, every 4-tuple exactly once  OK\n");
        else { say("    (A) FAILED at %u %u %u %u\n", failpos[0], failpos[1], failpos[2], failpos[3]); ok = false; }

        // (B) sharpness and (C) masking: all sets in full mode, every 40th / 20th otherwise
        unsigned pos5[K + 1] = {0, 1, 2, 3, 4};
        unsigned n5 = 0, idx = 0;
        do {
            if (full || idx++ % 40 == 0) {
                ++n5;
                if (count_reachable<F, M>(pos5, K + 1) != (1ull << (M * K))) {
                    say("    (B) FAILED at %u %u %u %u %u\n", pos5[0], pos5[1], pos5[2], pos5[3], pos5[4]);
                    ok = false; break;
                }
            }
        } while (next_combination(pos5, K + 1, F::SIZE));
        say("    (B) sharpness: %u/4368 sets of 5 positions reach 2^%u of 2^%u tuples  OK\n", n5, M * K, M * (K + 1));

        unsigned posm[K] = {0, 1, 2, 3};
        unsigned nm = 0; idx = 0;
        do {
            if (full || idx++ % 20 == 0) {
                ++nm;
                if (!check_exactly_once<F, M>(posm, K, true)) {
                    say("    (C) FAILED (masked) at %u %u %u %u\n", posm[0], posm[1], posm[2], posm[3]);
                    ok = false; break;
                }
            }
        } while (next_combination(posm, K, F::SIZE));
        say("    (C) masking:   XOR with fixed mask preserves (A) on %u/1820 sets  OK\n", nm);

        // (D) negative control: reducible x^4+1 must break (A)
        const bool control_fails = !gf16_check_A<0x1>(nullptr) && !SmallField<4, 0x1>::is_field();
        say("    (D) control:   reducible x^4+1 breaks (A)  %s\n", control_fails ? "OK" : "FAILED (test is not discriminating!)");
        ok &= control_fails;
    }

    if (full) {
        using F = SmallField<8, 0x1B>;
        constexpr unsigned M = 8, K = 3;
        say("  GF(2^8) k=%u seed space 2^%u%s\n", K, M * K, F::is_field() ? "" : "  FIELD CHECK FAILED");
        ok &= F::is_field();

        uint64_t rng = 0xD1B54A32D192ED03ULL;
        auto next = [&rng]() { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; };
        std::vector<std::vector<unsigned>> sets = {{0, 1, 2}, {0, 1, 255}, {0, 128, 255}, {1, 2, 4}};
        while (sets.size() < gf8_sets) {
            unsigned p[K];
            do { for (unsigned& v : p) v = unsigned(next() & 0xFF); }
            while (p[0] == p[1] || p[0] == p[2] || p[1] == p[2]);
            sets.push_back({p[0], p[1], p[2]});
        }
        if (sets.size() > gf8_sets) sets.resize(gf8_sets ? gf8_sets : 1);

        for (const auto& s : sets)
            if (!check_exactly_once<F, M>(s.data(), K, false)) {
                say("    (A) FAILED at %u %u %u\n", s[0], s[1], s[2]); ok = false; break;
            }
        say("    (A) k-wise:    %zu position sets, every 3-tuple exactly once  OK\n", sets.size());

        const unsigned pos4[K + 1] = {0, 1, 2, 3};
        if (count_reachable<F, M>(pos4, K + 1) != (1ull << (M * K))) { say("    (B) FAILED\n"); ok = false; }
        else say("    (B) sharpness: 4 positions reach 2^%u of 2^%u tuples  OK\n", M * K, M * (K + 1));

        if (!check_exactly_once<F, M>(sets[0].data(), K, true)) { say("    (C) FAILED (masked)\n"); ok = false; }
        else say("    (C) masking:   XOR with fixed mask preserves (A)  OK\n");
    }

    say("  toy test time: %.1f s\n", seconds_since(t0));
    return ok;
}
