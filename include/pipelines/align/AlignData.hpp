#pragma once

// Standard
#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

// Internal
#include "AlignSample.hpp"
#include "PreprocessData.hpp"
#include "pipelines/PipelineData.hpp"

namespace pipelines::align {

static const std::array<std::string, 1> validInputSuffixSingleton = {
    preprocess::outSampleFastqSuffix};
static const std::array<std::string, 5> validInputSuffixesPaired{
    preprocess::outSampleFastqPairedMergeSuffix,
    preprocess::outSampleFastqPairedForwardSingletonSuffix,
    preprocess::outSampleFastqPairedReverseSingletonSuffix,
    preprocess::outSampleFastqPairedForwardPairedSuffix,
    preprocess::outSampleFastqPairedReversePairedSuffix};

static const std::array<std::string, 6> validInputSuffixes{
    preprocess::outSampleFastqSuffix,
    preprocess::outSampleFastqPairedMergeSuffix,
    preprocess::outSampleFastqPairedForwardSingletonSuffix,
    preprocess::outSampleFastqPairedReverseSingletonSuffix,
    preprocess::outSampleFastqPairedForwardPairedSuffix,
    preprocess::outSampleFastqPairedReversePairedSuffix};

static const std::string outSampleAlignedSuffix = "_complete_aligned.bam";
static const std::string outSampleMergedAlignedSuffix = "_merged_aligned.bam";
static const std::string outSampleSingletonForwardAlignedSuffix = "_singleton_forward_aligned.bam";
static const std::string outSampleSingletonReverseAlignedSuffix = "_singleton_reverse_aligned.bam";
static const std::string outSamplePairedAlignedSuffix = "_paired_aligned.bam";

static const std::string pipelinePrefix = "02_align";

struct AlignData : public pipelines::PipelineData {
    std::vector<AlignSampleType> treatmentSamples;
    std::optional<std::vector<AlignSampleType>> controlSamples;

    AlignData(const fs::path& outputDir, const fs::path& treatmentDir,
              const std::optional<fs::path> controlDir)
        : treatmentSamples(retrieveSamples(treatmentSampleGroup, treatmentDir, outputDir)),
          controlSamples(controlDir ? std::optional(retrieveSamples(controlSampleGroup,
                                                                    controlDir.value(), outputDir))
                                    : std::nullopt) {};

   private:
    static auto retrieveSamples(const std::string& sampleGroup, const fs::path& parentDir,
                                const fs::path& outputDir) -> std::vector<AlignSampleType>;
    static auto retrieveInputSamples(const fs::path& parentDir) -> std::vector<InputSampleType>;
    static auto retrieveInputSample(const fs::path& sampleDir) -> InputSampleType;
    static auto retrieveInputPaired(const std::string& sampleName,
                                    const std::vector<fs::path>& inputSamples) -> AlignInputPaired;
    static auto retrieveInputSingle(const std::string& sampleName, const fs::path& inputSample)
        -> AlignInputSingle;
};

}  // namespace pipelines::align
