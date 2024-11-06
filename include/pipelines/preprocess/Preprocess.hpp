#pragma once

// Standard
#include <cstddef>
#include <ranges>
#include <vector>

// boost
#include <boost/program_options.hpp>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>

// Class
#include "Adapter.hpp"
#include "PreprocessData.hpp"
#include "PreprocessParameters.hpp"
#include "PreprocessSample.hpp"

using seqan3::operator""_dna5;

namespace pipelines::preprocess {
class Preprocess {
   public:
    explicit Preprocess(PreprocessParameters params);

    void process(const PreprocessData &data) const;

   private:
    using SingleEndAsyncInputBuffer = seqan3::detail::async_input_buffer_view<
        std::ranges::ref_view<seqan3::sequence_file_input<>>>;
    using PairedEndAsyncInputBuffer = seqan3::detail::async_input_buffer_view<std::views::all_t<
        seqan::stl::ranges::zip_view<std::ranges::ref_view<seqan3::sequence_file_input<>>,
                                     std::ranges::ref_view<seqan3::sequence_file_input<>>>>>;

    PreprocessParameters parameters;

    struct SingleEndResult {
        size_t passedRecords{0};
        size_t failedRecords{0};

        void operator+=(const SingleEndResult &other) {
            passedRecords += other.passedRecords;
            failedRecords += other.failedRecords;
        }
    };

    struct PairedEndResult {
        size_t mergedRecords{0};
        size_t singleFwdRecords{0};
        size_t singleRevRecords{0};
        size_t pairedRecordPairs{0};
        size_t failedMergedRecords{0};
        size_t failedForwardRecords{0};
        size_t failedReverseRecords{0};

        void operator+=(const PairedEndResult &other) {
            mergedRecords += other.mergedRecords;
            singleFwdRecords += other.singleFwdRecords;
            singleRevRecords += other.singleRevRecords;
            pairedRecordPairs += other.pairedRecordPairs;
            failedMergedRecords += other.failedMergedRecords;
            failedForwardRecords += other.failedForwardRecords;
            failedReverseRecords += other.failedReverseRecords;
        }
    };

    auto passesFilters(const auto &record) const -> bool;

    void processSample(const PreprocessSampleType &sample) const;

    void processSingleEnd(const PreprocessSampleSingle &sample) const;
    void processPairedEnd(const PreprocessSamplePaired &sample) const;

    auto processSingleEndRecordChunk(SingleEndAsyncInputBuffer &asyncInputBuffer,
                                     const std::vector<Adapter> &adapters5,
                                     const std::vector<Adapter> &adapters3,
                                     const fs::path &tmpOutDir) const -> SingleEndResult;

    auto processPairedEndRecordChunk(Preprocess::PairedEndAsyncInputBuffer &pairedRecordInputBuffer,
                                     const std::vector<Adapter> &adapters5f,
                                     const std::vector<Adapter> &adapters3f,
                                     const std::vector<Adapter> &adapters5r,
                                     const std::vector<Adapter> &adapters3r,
                                     const PrepocessSampleOutputPaired &sampleOutput) const
        -> PairedEndResult;
};

}  // namespace pipelines::preprocess
