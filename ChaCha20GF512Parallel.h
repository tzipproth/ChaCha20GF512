#pragma once

#include "ChaCha20GF512FFT.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <thread>
#include <vector>

// Bulk parallel convenience wrapper for ChaCha20GF512FFT.
//
// It does NOT create multiple RNG streams and does NOT combine streams.
// Instead it evaluates disjoint position ranges of the SAME logical stream in
// parallel. fill_parallel_at(..., start, ...) is bit-identical to a single
// ChaCha20GF512FFT seek(start) followed by the same number of next_int() calls.
//
// The middle of the requested range is split only at 1024-word GF/FFT output
// boundaries, avoiding redundant additive-FFT blocks between workers. A short unaligned
// prefix/suffix is handled serially.
class ChaCha20GF512Parallel
{
public:
    using Position = ChaCha20GF512FFT::Position;
    static constexpr std::size_t FFT_SIZE = ChaCha20GF512FFT::FFT_SIZE;
    static constexpr std::size_t OUTPUT_WORDS_PER_FFT = ChaCha20GF512FFT::OUTPUT_WORDS_PER_FFT;

    ChaCha20GF512Parallel(const std::uint8_t* chacha_seed32,
                         const std::uint8_t* gf_seed8192,
                         std::uint64_t stream_id = 0)
        : prototype_(chacha_seed32, gf_seed8192, stream_id)
    {
    }

    explicit ChaCha20GF512Parallel(const std::uint8_t* full_seed8224,
                                   std::uint64_t stream_id = 0)
        : prototype_(full_seed8224, stream_id)
    {
    }

    explicit ChaCha20GF512Parallel(std::uint64_t seed,
                                   std::uint64_t stream_id = 0)
        : prototype_(seed, stream_id)
    {
    }

    void seek(Position word_position)
    {
        position_ = word_position;
        exhausted_ = false;
    }

    void seek(std::uint64_t word_position)
    {
        seek(Position(word_position));
    }

    Position position128() const
    {
        return position_;
    }

    // Backward-compatible helper for callers that only use the old 64-bit range.
    std::uint64_t position() const
    {
        if (position_.hi != 0)
            throw std::overflow_error("ChaCha20GF512Parallel::position(): position exceeds 64 bits; use position128()");
        return position_.lo;
    }

    void fill_parallel(std::uint64_t* dst, std::size_t count,
                       unsigned thread_count)
    {
        if (count == 0)
            return;
        if (exhausted_)
            throw std::overflow_error("ChaCha20GF512Parallel: 128-bit word-position space exhausted");

        validate_range(position_, count);
        fill_parallel_at(dst, count, position_, thread_count);

        const std::uint64_t step = static_cast<std::uint64_t>(count);
        Position next = position_;
        if (next.add_u64(step)) {
            // The only valid overflow case is that the generated range ended
            // exactly at 2^128. There is then no representable next position.
            exhausted_ = true;
        }
        else {
            position_ = next;
        }
    }

    void fill_parallel_at(std::uint64_t* dst, std::size_t count,
                          Position start_position,
                          unsigned thread_count) const
    {
        if (count == 0)
            return;

        validate_range(start_position, count);

        if (thread_count == 0)
            thread_count = std::max(1u, std::thread::hardware_concurrency());

        // Small ranges are faster without thread creation overhead.
        if (thread_count <= 1 || count < OUTPUT_WORDS_PER_FFT * 2) {
            fill_serial_range(dst, count, start_position);
            return;
        }

        std::size_t written = 0;
        Position pos = start_position;

        // Reach a 1024-word output boundary without making two threads evaluate
        // the same 512-point GF FFT block.
        const std::size_t offset = pos.low10();
        if (offset != 0) {
            const std::size_t prefix = std::min(count, OUTPUT_WORDS_PER_FFT - offset);
            fill_serial_range(dst, prefix, pos);
            dst += prefix;
            written += prefix;
            add_checked(pos, prefix);
        }

        const std::size_t remaining = count - written;
        const std::size_t full_blocks = remaining / OUTPUT_WORDS_PER_FFT;
        const std::size_t middle_words = full_blocks * OUTPUT_WORDS_PER_FFT;

        if (full_blocks != 0) {
            const unsigned workers = static_cast<unsigned>(
                std::min<std::size_t>(thread_count, full_blocks));

            if (workers <= 1) {
                fill_serial_range(dst, middle_words, pos);
            }
            else {
                std::vector<std::thread> threads;
                threads.reserve(workers);
                std::vector<std::exception_ptr> errors(workers);

                const std::size_t base_blocks = full_blocks / workers;
                const std::size_t extra_blocks = full_blocks % workers;

                std::size_t block_cursor = 0;
                for (unsigned t = 0; t < workers; ++t) {
                    const std::size_t blocks = base_blocks + (t < extra_blocks ? 1 : 0);
                    const std::size_t word_offset = block_cursor * OUTPUT_WORDS_PER_FFT;
                    const std::size_t words = blocks * OUTPUT_WORDS_PER_FFT;
                    Position worker_pos = pos;
                    add_checked(worker_pos, word_offset);
                    std::uint64_t* worker_dst = dst + word_offset;

                    threads.emplace_back([this, worker_dst, words, worker_pos, &errors, t]() {
                        try {
                            fill_serial_range(worker_dst, words, worker_pos);
                        }
                        catch (...) {
                            errors[t] = std::current_exception();
                        }
                    });
                    block_cursor += blocks;
                }

                for (auto& th : threads)
                    th.join();
                for (const auto& error : errors)
                    if (error) std::rethrow_exception(error);
            }

            dst += middle_words;
            written += middle_words;
            add_checked(pos, middle_words);
        }

        const std::size_t suffix = count - written;
        if (suffix != 0)
            fill_serial_range(dst, suffix, pos);
    }

    void fill_parallel_at(std::uint64_t* dst, std::size_t count,
                          std::uint64_t start_position,
                          unsigned thread_count) const
    {
        fill_parallel_at(dst, count, Position(start_position), thread_count);
    }

private:
    ChaCha20GF512FFT prototype_;
    Position position_{};
    bool exhausted_ = false;

    static void add_checked(Position& p, std::size_t amount)
    {
        if (amount == 0)
            return;
        if (p.add_u64(static_cast<std::uint64_t>(amount)))
            throw std::overflow_error("ChaCha20GF512Parallel: position overflow");
    }

    static void validate_range(Position start, std::size_t count)
    {
        // Validate the last requested word, not the exclusive end. This allows
        // a one-word request at position 2^128-1.
        if (count == 0)
            return;
        Position last = start;
        const std::uint64_t delta = static_cast<std::uint64_t>(count - 1);
        if (last.add_u64(delta))
            throw std::overflow_error("ChaCha20GF512Parallel: requested range exceeds 128-bit position space");
    }

    void fill_serial_range(std::uint64_t* dst, std::size_t count,
                           Position start_position) const
    {
        // Copying roughly 20-30 KiB once per bulk worker is cheap relative to
        // processing a large range and avoids any shared mutable RNG state.
        ChaCha20GF512FFT worker = prototype_;
        worker.seek(start_position);
        worker.generate(dst, count);
    }
};
