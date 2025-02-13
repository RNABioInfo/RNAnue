#pragma once

// Standard
#include <cstddef>
#include <deque>
#include <filesystem>
#include <memory>
#include <optional>
#include <ranges>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/io/sam_file/input.hpp>

// Internal
#include "AsyncSplitRecordGroupBuffer.hpp"
#include "DetectData.hpp"
#include "DetectParameters.hpp"
#include "DetectSample.hpp"
#include "FeatureAnnotator.hpp"
#include "SamRecord.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"
#include "SplitRecordsEvaluator.hpp"

using namespace dataTypes;

namespace pipelines::detect {

using namespace annotation;
namespace fs = std::filesystem;

class Detect {
   public:
    explicit Detect(DetectParameters params) : params(std::move(params)) {};

    Detect(const Detect &) = default;
    Detect(Detect &&) = delete;
    auto operator=(const Detect &) -> Detect & = delete;
    auto operator=(Detect &&) -> Detect & = delete;
    ~Detect() = default;

    void process(const DetectData &data);

   private:
    using AsyncGroupBufferType = AsyncSplitRecordGroupBufferView<std::ranges::ref_view<
        seqan3::sam_file_input<seqan3::sam_file_input_default_traits<>, SamFieldIDs>>>;

    struct ChunkedOutTmpDirs {
        fs::path outputTmpSplitsDir;
        fs::path outputTmpMultisplitsDir;
        fs::path outputTmpUnassignedContiguousDir;
    };

    using TranscriptCounts = std::unordered_map<std::string, size_t>;

    struct Result {
        size_t processedRecordsCount{0};
        TranscriptCounts transcriptCounts;
        size_t splitFragmentsCount{0};
        size_t singletonFragmentsCount{0};
        size_t removedDueToLowMappingQuality{0};
        size_t removedDueToFragmentLength{0};

        void operator+=(const Result &other) {
            processedRecordsCount += other.processedRecordsCount;
            splitFragmentsCount += other.splitFragmentsCount;
            singletonFragmentsCount += other.singletonFragmentsCount;
            removedDueToLowMappingQuality += other.removedDueToLowMappingQuality;
            removedDueToFragmentLength += other.removedDueToFragmentLength;

            for (const auto &[transcript, count] : other.transcriptCounts) {
                transcriptCounts[transcript] += count;
            }
        }
    };

    DetectParameters params;

    std::optional<SplitRecordsEvaluator> splitRecordsEvaluator;

    [[nodiscard]] static auto getSplitRecordsEvaluatorParameters(
        const DetectParameters &params, std::shared_ptr<const FeatureAnnotator> featureAnnotator)
        -> SplitRecordsEvaluationParameters::ParameterVariant;

    static auto getReferenceIDs(const fs::path &mappingsInPath) -> std::deque<std::string>;

    void processSample(const DetectSample &sample, std::shared_ptr<const FeatureAnnotator>) const;

    auto processRecordChunk(const ChunkedOutTmpDirs &outTmpDirs,
                            AsyncGroupBufferType &recordInputBuffer,
                            const std::deque<std::string> &refIDs,
                            const std::vector<size_t> &refLengths,
                            std::shared_ptr<const FeatureAnnotator> featureAnnotator) const
        -> Result;
    auto processReadRecords(const std::vector<SamRecord> &readRecords, auto &splitsOut,
                            [[maybe_unused]] auto &multiSplitsOut) const -> size_t;

    [[nodiscard]] auto constructSplitRecords(const SamRecord &readRecord) const
        -> std::optional<SplitRecords>;
    [[nodiscard]] auto constructSplitRecords(const std::vector<SamRecord> &readRecords) const
        -> std::optional<SplitRecords>;
    [[nodiscard]] auto getBestSplitRecords(const std::vector<SamRecord> &readRecords) const
        -> std::optional<SplitRecordsEvaluator::EvaluatedSplitRecords>;

    static void mergeOutputFiles(const ChunkedOutTmpDirs &tmpDirs, const DetectOutput &output);
    static void writeSamFile(auto &samOut, const std::vector<SamRecord> &splitRecords);

    static void writeTranscriptCountsFile(const fs::path &transcriptCountsFilePath,
                                          const TranscriptCounts &transcriptCounts);
    static void writeReadCountsSummaryFile(const Result &results, const std::string &sampleName,
                                           const fs::path &statsFilePath);

    [[nodiscard]] static auto prepareTmpOutputDirs(const fs::path &tmpOutDir) -> ChunkedOutTmpDirs;
};

}  // namespace pipelines::detect
