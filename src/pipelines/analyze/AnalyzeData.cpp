#include "AnalyzeData.hpp"

// Standard
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

// Internal
#include "AnalyzeSample.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "Utility.hpp"

using namespace helper;

namespace pipelines::analyze {

[[nodiscard]] auto AnalyzeData::retrieveSamples(const std::string& sampleGroup,
                                                const fs::path& parentDir,
                                                const fs::path& outputDir)
    -> std::vector<AnalyzeSample> {
    const std::vector<AnalyzeInput> inputSamples = retrieveInputSamples(parentDir);

    std::vector<AnalyzeSample> samples;
    samples.reserve(inputSamples.size());

    const fs::path outputDirPipeline = outputDir / pipelinePrefix / sampleGroup;

    if (fs::exists(outputDirPipeline)) {
        Logger::log<LogLevel::WARNING>(
            "Output directory already exists. Results will be overwritten.");
    }

    for (const AnalyzeInput& inputSample : inputSamples) {
        const fs::path outputDirSample = outputDirPipeline / inputSample.sampleName;

        fs::create_directories(outputDirSample);

        const fs::path interactionsPath =
            outputDirSample / (inputSample.sampleName + outInteractionsSuffix);
        const fs::path interactionReadIDsPath =
            outputDirSample / (inputSample.sampleName + outInteractionReadIDsSuffix);
        const fs::path interactionsTranscriptCountsPath =
            outputDirSample / (inputSample.sampleName + outInteractionsTranscriptCountsSuffix);
        const fs::path interactionsBEDPath =
            outputDirSample / (inputSample.sampleName + outInteractionsBEDSuffix);
        const fs::path interactionsBEDARCPath =
            outputDirSample / (inputSample.sampleName + outInteractionsBEDARCSuffix);
        const fs::path supplementaryFeaturesPath =
            outputDirSample / (inputSample.sampleName + outSupplementaryFeaturesSuffix);

        samples.push_back(AnalyzeSample(
            inputSample,
            AnalyzeOutput{.interactionsPath = interactionsPath,
                          .interactionsReadIDsPath = interactionReadIDsPath,
                          .interactionsTranscriptCountsPath = interactionsTranscriptCountsPath,
                          .interactionsBEDPath = interactionsBEDPath,
                          .interactionsBEDARCPath = interactionsBEDARCPath,
                          .supplementaryFeaturesPath = supplementaryFeaturesPath}));
    }

    return samples;
}

[[nodiscard]] auto AnalyzeData::retrieveInputSamples(const fs::path& parentDir)
    -> std::vector<AnalyzeInput> {
    const std::vector<fs::path> sampleDirs = getSubDirectories(parentDir);

    std::vector<AnalyzeInput> samples;

    for (const fs::path& sampleDir : sampleDirs) {
        const std::vector<fs::path> sampleFiles = getValidFilePaths(sampleDir, validSuffices);

        if (sampleFiles.size() != validSuffices.size()) {
            const std::string message = "Expected " + std::to_string(validSuffices.size()) +
                                        " files in " + sampleDir.string() + " but found " +
                                        std::to_string(sampleFiles.size()) + " valid files";
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(message);
            throw std::runtime_error(message);
        }

        const std::string sampleName = getSampleName(sampleFiles.front());

        std::optional<fs::path> splitAlignmentsPath;
        std::optional<fs::path> multisplitAlignmentsPath;
        std::optional<fs::path> unassignedContiguousAlignmentsPath;
        std::optional<fs::path> contiguousAlignmentsTranscriptCountsPath;
        std::optional<fs::path> sharedReadCountsPath;

        for (const fs::path& sampleFile : sampleFiles) {
            if (hasSuffix(sampleFile, validInputSplitAlignmentsSuffix)) {
                splitAlignmentsPath = sampleFile;
            } else if (hasSuffix(sampleFile, validInputMultisplitAlignmentsSuffix)) {
                multisplitAlignmentsPath = sampleFile;
            } else if (hasSuffix(sampleFile, validInputUnassignedContiguousAlignmentsSuffix)) {
                unassignedContiguousAlignmentsPath = sampleFile;
            } else if (hasSuffix(sampleFile,
                                 validInputContiguousAlignmentsTranscriptCountsSuffix)) {
                contiguousAlignmentsTranscriptCountsPath = sampleFile;
            } else if (hasSuffix(sampleFile, validSharedReadCountsSuffix)) {
                sharedReadCountsPath = sampleFile;
            } else {
                const std::string message =
                    "Unexpected file " + sampleFile.string() + " found in " + sampleDir.string();
                Logger::log<LogLevel::WARNING>(message);
            }
        }

        if (!splitAlignmentsPath.has_value() || !multisplitAlignmentsPath.has_value() ||
            !unassignedContiguousAlignmentsPath.has_value() ||
            !contiguousAlignmentsTranscriptCountsPath.has_value() ||
            !sharedReadCountsPath.has_value()) {
            const std::string message =
                "Missing one or more required files in " + sampleDir.string();
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(message);
            throw std::runtime_error(message);
        }

        samples.push_back(AnalyzeInput{
            .sampleName = sampleName,
            .splitAlignmentsPath = splitAlignmentsPath.value(),
            .multisplitAlignmentsPath = multisplitAlignmentsPath.value(),
            .unassignedContiguousAlignmentsPath = unassignedContiguousAlignmentsPath.value(),
            .contiguousAlignmentsTranscriptCountsPath =
                contiguousAlignmentsTranscriptCountsPath.value(),
            .sampleFragmentCountsPath = sharedReadCountsPath.value()});
    }

    return samples;
}

auto AnalyzeData::getInputFilePaths() const -> std::vector<fs::path> {
    std::vector<fs::path> inputFilePaths;
    inputFilePaths.reserve(treatmentSamples.size());

    for (const auto& sample : treatmentSamples) {
        inputFilePaths.push_back(sample.input.splitAlignmentsPath);
    }

    if (controlSamples) {
        for (const auto& sample : controlSamples.value()) {
            inputFilePaths.push_back(sample.input.splitAlignmentsPath);
        }
    }

    return inputFilePaths;
}
}  // namespace pipelines::analyze
