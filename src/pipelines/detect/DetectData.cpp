#include "DetectData.hpp"

#include "Logger.hpp"
#include "Utility.hpp"

using namespace helper;

namespace pipelines::detect {

auto DetectData::retrieveSamples(const std::string& sampleGroup, const fs::path& parentDir,
                                 const fs::path& outputDir) -> std::vector<DetectSample> {
    const std::vector<DetectInput> inputSamples = retrieveInputSamples(parentDir);

    std::vector<DetectSample> samples;
    samples.reserve(inputSamples.size());

    const fs::path outputDirPipeline = outputDir / pipelinePrefix / sampleGroup;

    for (const DetectInput& inputSample : inputSamples) {
        const fs::path outputDirSample = outputDirPipeline / inputSample.sampleName;

        fs::create_directories(outputDirSample);

        const fs::path outputSplitAlignmentsPath =
            outputDirSample / (inputSample.sampleName + outSampleSplitAlignmentsSuffix);
        const fs::path outputMultisplitAlignmentsPath =
            outputDirSample / (inputSample.sampleName + outSampleMultisplitAlignmentsSuffix);
        const fs::path outputUnassignedContiguousAlignmentsPath =
            outputDirSample /
            (inputSample.sampleName + outSampleUnassignedContiguousAlignmentsSuffix);
        const fs::path outputContiguousAlignmentsTranscriptCountsPath =
            outputDirSample /
            (inputSample.sampleName + outSampleContiguousAlignmentsTranscriptCountsSuffix);

        const fs::path outputSharedStatsPath =
            outputDirSample / (inputSample.sampleName + outSampleReadCountsSummarySuffix);

        samples.push_back(DetectSample(
            inputSample,
            DetectOutput{outputSplitAlignmentsPath, outputMultisplitAlignmentsPath,
                         outputUnassignedContiguousAlignmentsPath,
                         outputContiguousAlignmentsTranscriptCountsPath, outputSharedStatsPath}));
    }

    return samples;
}

auto DetectData::retrieveInputSamples(const fs::path& parentDir) -> std::vector<DetectInput> {
    const std::vector<fs::path> sampleDirs = getSubDirectories(parentDir);

    std::vector<DetectInput> samples;

    for (const fs::path& sampleDir : sampleDirs) {
        const std::vector<fs::path> sampleFiles = getValidFilePaths(sampleDir, {validInputSuffix});

        if (sampleFiles.size() != 1) {
            const std::string message = "Expected 1 input file in " + sampleDir.string() +
                                        " but found " + std::to_string(sampleFiles.size()) +
                                        " valid files";
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(message);
            throw std::runtime_error(message);
        }

        const fs::path& sampleFile = sampleFiles.front();
        const auto sampleName = getSampleName(sampleFile);

        if (!hasSuffix(sampleFile, validInputSuffix)) {
            const std::string message = "Invalid input file suffix: " + sampleFile.string() +
                                        " . Expected: " + validInputSuffix;
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(message);
            throw std::runtime_error(message);
        }

        Logger::log<LogLevel::DEBUG>("Found input sample: " + sampleName +
                                     " file: " + sampleFile.string());

        samples.emplace_back(sampleName, sampleFile);
    }

    return samples;
}

}  // namespace pipelines::detect
