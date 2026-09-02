#pragma once

#include "ChaCha20GF512FFT.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

#if defined(_WIN32)
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <windows.h>
#  include <bcrypt.h>
#  if defined(_MSC_VER)
#    pragma comment(lib, "bcrypt.lib")
#  endif
#elif defined(__linux__)
#  include <cerrno>
#  include <sys/random.h>
#else
#  error "ChaCha20GF512Seed.h currently supports Windows and Linux."
#endif

// Convenience helpers for practical OS-backed initialization.
//
// IMPORTANT:
// BCryptGenRandom/getrandom() are operating-system CSPRNG interfaces. They are
// the recommended practical way to seed the generator, but reading 4128 bytes
// from an OS CSPRNG does NOT prove that 33024 bits of fresh physical entropy
// were injected. The exact information-theoretic 512-wise statement is a
// statement about uniform sampling of the full explicit seed space.
namespace chacha20gf512_seed
{
    inline bool os_random_bytes(void* dst, std::size_t bytes)
    {
        auto* p = static_cast<std::uint8_t*>(dst);

#if defined(_WIN32)
        // BCryptGenRandom takes ULONG lengths. Chunk for completeness.
        while (bytes != 0) {
            const ULONG chunk = static_cast<ULONG>(
                bytes > static_cast<std::size_t>(0xFFFFFFFFu)
                    ? 0xFFFFFFFFu : bytes);
            const NTSTATUS st = BCryptGenRandom(
                nullptr, reinterpret_cast<PUCHAR>(p), chunk,
                BCRYPT_USE_SYSTEM_PREFERRED_RNG);
            if (st < 0)
                return false;
            p += chunk;
            bytes -= chunk;
        }
        return true;

#elif defined(__linux__)
        while (bytes != 0) {
            const ssize_t n = ::getrandom(p, bytes, 0);
            if (n > 0) {
                p += static_cast<std::size_t>(n);
                bytes -= static_cast<std::size_t>(n);
                continue;
            }
            if (n < 0 && errno == EINTR)
                continue;
            return false;
        }
        return true;
#endif
    }

    using FullSeed = std::array<std::uint8_t, ChaCha20GF512FFT::FULL_SEED_BYTES>;

    inline FullSeed make_os_full_seed()
    {
        FullSeed seed{};
        if (!os_random_bytes(seed.data(), seed.size()))
            throw std::runtime_error("OS random generator failed");
        return seed;
    }

    inline ChaCha20GF512FFT make_os_seeded_rng(std::uint64_t stream_id = 0)
    {
        const FullSeed seed = make_os_full_seed();
        return ChaCha20GF512FFT(seed.data(), stream_id);
    }
}
