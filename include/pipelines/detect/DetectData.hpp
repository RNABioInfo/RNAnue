#pragma once

// Standard
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

// Internal
#include "AlignData.hpp"
#include "DetectSample.hpp"
#include "PipelineData.hpp"

namespace pipelines::detect {
namespace fs = std::filesystem;

static const std::string validInputSuffix = pipelines::align::outSampleAlignedSuffix;

static const std::string outSampleSplitAlignmentsSuffix = "_splits.bam";
static const std::string outSampleMultisplitAlignmentsSuffix = "_multisplits.bam";
static const std::string outSampleUnassignedContiguousAlignmentsSuffix =
    "_contiguous_unassigned_aligned.bam";
static const std::string outSampleContiguousAlignmentsTranscriptCountsSuffix =
    "_contiguous_transcript_counts.tsv";
static const std::string outSampleReadCountsSummarySuffix = "_read_counts_summary.tsv";

static const std::string pipelinePrefix = "03_detect";

struct DetectData : public PipelineData {
    std::vector<DetectSample> treatmentSamples;
    std::optional<std::vector<DetectSample>> controlSamples;

    DetectData(const fs::path& outputDir, const fs::path& treatmentDir,
               const std::optional<fs::path>& controlDir)
        : treatmentSamples(retrieveSamples(treatmentSampleGroup, treatmentDir, outputDir)),
          controlSamples(controlDir ? std::optional(retrieveSamples(controlSampleGroup,
                                                                    controlDir.value(), outputDir))
                                    : std::nullopt) {};

    [[nodiscard]] auto getInputFilePaths() const -> std::vector<fs::path>;

   private:
    static auto retrieveSamples(const std::string& sampleGroup, const fs::path& parentDir,
                                const fs::path& outputDir) -> std::vector<DetectSample>;

    static auto retrieveInputSamples(const fs::path& parentDir) -> std::vector<DetectInput>;
};

}  // namespace pipelines::detect
