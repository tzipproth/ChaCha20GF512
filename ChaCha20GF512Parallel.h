#pragma once

#include "ChaCha20GF512FFT.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>

// Bulk parallel convenience wrapper for ChaCha20GF512FFT.
//
// It does NOT create multiple RNG streams and does NOT combine streams.
// Instead it evaluates disjoint position ranges of the SAME logical stream in
// parallel. Therefore fill_parallel_at(..., start, ...) is bit-identical to a
// single ChaCha20GF512FFT seek(start) followed by the same number of next_int()
// calls.
//
// The middle of the requested range is split only at 512-word FFT boundaries,
// avoiding redundant additive-FFT blocks between workers. A short unaligned
// prefix/suffix is handled serially.
class ChaCha20GF512Parallel
{
public:
    static constexpr std::size_t FFT_SIZE = ChaCha20GF512FFT::FFT_SIZE;

    ChaCha20GF512Parallel(const std::uint8_t* chacha_seed32,
                        const std::uint8_t* gf_seed4096,
                        std::uint64_t stream_id = 0)
        : prototype_(chacha_seed32, gf_seed4096, stream_id)
    {
    }

    explicit ChaCha20GF512Parallel(const std::uint8_t* full_seed4128,
                                 std::uint64_t stream_id = 0)
        : prototype_(full_seed4128, stream_id)
    {
    }

    explicit ChaCha20GF512Parallel(std::uint64_t seed,
                                 std::uint64_t stream_id = 0)
        : prototype_(seed, stream_id)
    {
    }

    void seek(std::uint64_t word_position)
    {
        position_ = word_position;
    }

    std::uint64_t position() const
    {
        return position_;
    }

    void fill_parallel(std::uint64_t* dst, std::size_t count,
                       unsigned thread_count)
    {
        fill_parallel_at(dst, count, position_, thread_count);
        position_ += static_cast<std::uint64_t>(count);
    }

    void fill_parallel_at(std::uint64_t* dst, std::size_t count,
                          std::uint64_t start_position,
                          unsigned thread_count) const
    {
        if (count == 0)
            return;

        if (thread_count == 0)
            thread_count = std::max(1u, std::thread::hardware_concurrency());

        // Small ranges are faster without thread creation overhead.
        if (thread_count <= 1 || count < FFT_SIZE * 2) {
            fill_serial_range(dst, count, start_position);
            return;
        }

        std::size_t written = 0;
        std::uint64_t pos = start_position;

        // Reach a 512-word boundary without making two threads evaluate the
        // same FFT block.
        const std::size_t offset = static_cast<std::size_t>(pos & (FFT_SIZE - 1));
        if (offset != 0) {
            const std::size_t prefix = std::min(count, FFT_SIZE - offset);
            fill_serial_range(dst, prefix, pos);
            dst += prefix;
            written += prefix;
            pos += static_cast<std::uint64_t>(prefix);
        }

        const std::size_t remaining = count - written;
        const std::size_t full_blocks = remaining / FFT_SIZE;
        const std::size_t middle_words = full_blocks * FFT_SIZE;

        if (full_blocks != 0) {
            const unsigned workers = static_cast<unsigned>(
                std::min<std::size_t>(thread_count, full_blocks));

            if (workers <= 1) {
                fill_serial_range(dst, middle_words, pos);
            }
            else {
                std::vector<std::thread> threads;
                threads.reserve(workers);

                const std::size_t base_blocks = full_blocks / workers;
                const std::size_t extra_blocks = full_blocks % workers;

                std::size_t block_cursor = 0;
                for (unsigned t = 0; t < workers; ++t) {
                    const std::size_t blocks = base_blocks + (t < extra_blocks ? 1 : 0);
                    const std::size_t word_offset = block_cursor * FFT_SIZE;
                    const std::size_t words = blocks * FFT_SIZE;
                    const std::uint64_t worker_pos =
                        pos + static_cast<std::uint64_t>(word_offset);
                    std::uint64_t* worker_dst = dst + word_offset;

                    threads.emplace_back([this, worker_dst, words, worker_pos]() {
                        fill_serial_range(worker_dst, words, worker_pos);
                    });
                    block_cursor += blocks;
                }

                for (auto& th : threads)
                    th.join();
            }

            dst += middle_words;
            written += middle_words;
            pos += static_cast<std::uint64_t>(middle_words);
        }

        const std::size_t suffix = count - written;
        if (suffix != 0)
            fill_serial_range(dst, suffix, pos);
    }

private:
    ChaCha20GF512FFT prototype_;
    std::uint64_t position_ = 0;

    void fill_serial_range(std::uint64_t* dst, std::size_t count,
                           std::uint64_t start_position) const
    {
        // Copying ~10 KB once per bulk worker is cheap relative to processing a
        // large range and avoids any shared mutable RNG state.
        ChaCha20GF512FFT worker = prototype_;
        worker.seek(start_position);
        for (std::size_t i = 0; i < count; ++i)
            dst[i] = worker.next_int();
    }
};
