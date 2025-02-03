#pragma once

// standard
#include <cstdlib>
#include <ranges>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/io/views/async_input_buffer.hpp>

// Internal
#include "Adapter.hpp"
#include "FastqRecord.hpp"
#include "PreprocessParameters.hpp"
#include "PreprocessSample.hpp"
#include "RecordTrimmer.hpp"
#include "seqan3/contrib/std/zip_view.hpp"
#include "seqan3/io/sequence_file/input.hpp"

namespace pipelines::preprocess {

using namespace dataTypes;

class PairedEndPreprocessor {
   public:
    PairedEndPreprocessor(PreprocessParameters parameters) : parameters(std::move(parameters)) {}

    void process(const PreprocessSamplePaired& sample) const;

   private:
    using PairedEndAsyncInputBuffer = seqan3::detail::async_input_buffer_view<std::views::all_t<
        seqan::stl::ranges::zip_view<std::ranges::ref_view<seqan3::sequence_file_input<>>,
                                     std::ranges::ref_view<seqan3::sequence_file_input<>>>>>;

    struct ChunkResult {
        [[nodiscard]] auto getMergedRecords() const { return mergedRecords; }
        [[nodiscard]] auto getSingleFwdRecords() const { return singleFwdRecords; }
        [[nodiscard]] auto getSingleRevRecords() const { return singleRevRecords; }
        [[nodiscard]] auto getPairedRecords() const { return pairedRecords; }
        [[nodiscard]] auto getFailedMergedRecords() const { return failedMergedRecords; }
        [[nodiscard]] auto getFailedForwardRecords() const { return failedForwardRecords; }
        [[nodiscard]] auto getFailedReverseRecords() const { return failedReverseRecords; }

        void incrementMergedRecords() { ++mergedRecords; }
        void incrementSingleFwdRecords() { ++singleFwdRecords; }
        void incrementSingleRevRecords() { ++singleRevRecords; }
        void incrementPairedRecords() { ++pairedRecords; }
        void incrementFailedMergedRecords() { ++failedMergedRecords; }
        void incrementFailedForwardRecords() { ++failedForwardRecords; }
        void incrementFailedReverseRecords() { ++failedReverseRecords; }

        void operator+=(const PairedEndPreprocessor::ChunkResult& other);

       private:
        size_t mergedRecords{0};
        size_t singleFwdRecords{0};
        size_t singleRevRecords{0};
        size_t pairedRecords{0};
        size_t failedMergedRecords{0};
        size_t failedForwardRecords{0};
        size_t failedReverseRecords{0};
    };

    PreprocessParameters parameters;

    std::vector<Adapter> adapters5fwd;
    std::vector<Adapter> adapters3fwd;
    std::vector<Adapter> adapters5rev;
    std::vector<Adapter> adapters3rev;

    [[nodiscard]] auto processWithDeduplication(const PreprocessSamplePaired& sample) const
        -> ChunkResult;
    [[nodiscard]] auto processWithoutDeduplication(const PreprocessSamplePaired& sample) const
        -> ChunkResult;

    void trimWindowedQuality(PairedFastqRecords& records,
                             const RecordTrimmer::TrimWindowedConfig& config) const;

    void trimAdapters(PairedFastqRecords& records) const;

    template <typename T>
    auto processChunk(T& recordIterator, const PrepocessSampleOutputPaired& tmpOutDir) const
        -> ChunkResult;
};

}  // namespace pipelines::preprocess
