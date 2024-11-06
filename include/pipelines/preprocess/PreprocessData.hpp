#pragma once

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "PreprocessSample.hpp"
#include "pipelines/PipelineData.hpp"

namespace pipelines::preprocess {

namespace fs = std::filesystem;

static const std::array<std::string, 4> validInputSuffixes{".fastq", ".fastq.gz", ".fq", ".fq.gz"};
static const std::array<std::string, 3> validForwardTags{"_R1", "_1", "_forward"};
static const std::array<std::string, 3> validReverseTags{"_R2", "_2", "_reverse"};

static const std::string outSampleFastqSuffix = "_passed.fastq.gz";
static const std::string outSampleTmpFastqDirPrefix = "tmp_fastq";
static const std::string outSampleTmpMergedFastqDirPrefix = "tmp_merged_fastq";
static const std::string outSampleTmpForwardSingletonFastqDirPrefix = "tmp_forward_singleton_fastq";
static const std::string outSampleTmpReverseSingletonFastqDirPrefix = "tmp_reverse_singleton_fastq";
static const std::string outSampleTmpForwardPairedFastqDirPrefix = "tmp_forward_paired_fastq";
static const std::string outSampleTmpReversePairedFastqDirPrefix = "tmp_reverse_paired_fastq";

static const std::string outSampleFastqPairedMergeSuffix = "_merged_passed.fastq.gz";
static const std::string outSampleFastqPairedForwardSingletonSuffix =
    "_singleton_passed_R1.fastq.gz";
static const std::string outSampleFastqPairedReverseSingletonSuffix =
    "_singleton_passed_R2.fastq.gz";
static const std::string outSampleFastqPairedForwardPairedSuffix = "_non_merged_passed_R1.fastq.gz";
static const std::string outSampleFastqPairedReversePairedSuffix = "_non_merged_passed_R2.fastq.gz";

static const std::string pipelinePrefix = "01_preprocess";

struct PreprocessData : public pipelines::PipelineData {
    std::vector<PreprocessSampleType> treatmentSamples;
    std::optional<std::vector<PreprocessSampleType>> controlSamples;

    PreprocessData(const fs::path& outputDir, const fs::path& treatmentDir,
                   const std::optional<fs::path>& controlDir)
        : treatmentSamples(retrieveSamples(treatmentSampleGroup, treatmentDir, outputDir)),
          controlSamples(controlDir ? std::optional(retrieveSamples(controlSampleGroup,
                                                                    controlDir.value(), outputDir))
                                    : std::nullopt) {};

   private:
    static auto retrieveSamples(const std::string& sampleGroup, const fs::path& parentDir,
                                const fs::path& outputDir) -> std::vector<PreprocessSampleType>;
    static auto retrieveInputSamples(const fs::path& parentDir) -> std::vector<InputSampleType>;
    static auto retrieveInputSample(const fs::path& sampleDir) -> InputSampleType;
    static auto validatePairedFilePaths(std::pair<fs::path, fs::path>& pairedInputPaths) -> bool;
};

}  // namespace pipelines::preprocess
