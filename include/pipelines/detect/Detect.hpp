#pragma once

// Standard
#include <cstddef>
#include <deque>
#include <filesystem>
#include <format>
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
#include "AsyncSplitReadGroupBuffer.hpp"
#include "ComplementarityEvaluationStep.hpp"
#include "DetectData.hpp"
#include "DetectParameters.hpp"
#include "DetectSample.hpp"
#include "HybridizationEvaluationStep.hpp"
#include "ReadGroupEvaluationParameters.hpp"
#include "ReadGroupPostScoringStep.hpp"
#include "SamRecord.hpp"
#include "SamReference.hpp"
#include "SegemehlReadGroupPreprocessor.hpp"
#include "TempOutputDirs.hpp"

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
    using AsyncGroupBufferType = AsyncSplitReadGroupBufferView<std::ranges::ref_view<
        seqan3::sam_file_input<seqan3::sam_file_input_default_traits<>, SamFieldIDs>>>;

    using ReadGroup = std::vector<SamRecord>;

    using TranscriptCounts = std::unordered_map<std::string, double>;

    struct Result {
        SegemehlReadGroupPreprocessorMetrics preprocessMetrics;
        ComplementarityEvaluationStepMetrics complementarityMetrics;
        HybridizationEvaluationStepMetrics hybridizationMetrics;
        ReadGroupPostScoringStepMetrics postprocessMetrics;

        TranscriptCounts singletonTranscriptCounts;

        void createPlots(const fs::path &outDir) const noexcept {
            preprocessMetrics.createPlots(outDir);
            complementarityMetrics.createPlots(outDir);
            hybridizationMetrics.createPlots(outDir);
            postprocessMetrics.createPlots(outDir);
        }

        [[nodiscard]] auto toString() const -> std::string {
            return std::format(
                "\n\tPreprocessed {} read groups. Of which {} had at least one passed hit group "
                "and {} failed completely. \n\tComplementarity evaluation resulted in {} passed "
                "split hit groups and {} failed split hit groups. \n\tHybridization evalutaion "
                "resulted in {} passed split hit groups and {} failed split hit groups. "
                "\n\tPostprocessed {} read groups with a total contribution score of {}.",
                preprocessMetrics.getTotalReadGroupCount(),
                preprocessMetrics.getSuccesReadGroupCount(),
                preprocessMetrics.getTotalFailedReadGroupCount(),
                complementarityMetrics.getPassedCount(), complementarityMetrics.getFailedCount(),
                hybridizationMetrics.getPassedCount(), hybridizationMetrics.getFailedCount(),
                postprocessMetrics.getReadGroupCount(),
                postprocessMetrics.getContributionScoreSum());
        }

        void operator+=(const Result &other) {
            preprocessMetrics += other.preprocessMetrics;
            complementarityMetrics += other.complementarityMetrics;
            hybridizationMetrics += other.hybridizationMetrics;
            postprocessMetrics += other.postprocessMetrics;

            for (const auto &[transcript, count] : other.singletonTranscriptCounts) {
                singletonTranscriptCounts[transcript] += count;
            }
        }
    };

    DetectParameters params;

    static auto getReferenceIDs(const fs::path &mappingsInPath) -> std::deque<std::string>;

    template <ReadGroupEvaluationParameters::Type ParamT>
    void processSample(const DetectSample &sample, const ParamT &evaluationParams) const;

    template <ReadGroupEvaluationParameters::Type ParamT>
    auto processRecordChunk(const TempOutputDirs &outTmpDirs,
                            AsyncGroupBufferType &recordInputBuffer, SamReference &reference,
                            const ParamT &evaluationParams) const -> Result;

    static void mergeTmpFiles(const TempOutputDirs &tmpDirs, const DetectOutput &output,
                              const SamReference &reference);

    static void writeTranscriptCountsFile(const fs::path &transcriptCountsFilePath,
                                          const TranscriptCounts &transcriptCounts);
    static void writeReadCountsSummaryFile(const Result &results, const std::string &sampleName,
                                           const fs::path &statsFilePath);

    [[nodiscard]] static auto prepareTmpOutputDirs(const fs::path &tmpOutDir) -> TempOutputDirs;
};

}  // namespace pipelines::detect
