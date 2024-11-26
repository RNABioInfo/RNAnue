#include "AlignData.hpp"

#include <variant>
#include <vector>

#include "AlignSample.hpp"
#include "PreprocessData.hpp"
#include "Utility.hpp"

using namespace helper;

namespace pipelines::align {
auto AlignData::retrieveSamples(const std::string& sampleGroup, const fs::path& parentDir,
                                const fs::path& outputDir) -> std::vector<AlignSampleType> {
    const std::vector<InputSampleType> inputSamples = retrieveInputSamples(parentDir);

    std::vector<AlignSampleType> samples;
    samples.reserve(inputSamples.size());

    const fs::path outputDirPipeline = outputDir / pipelinePrefix / sampleGroup;

    for (const InputSampleType& inputSample : inputSamples) {
        if (const auto* inputSampleSingle = std::get_if<AlignInputSingle>(&inputSample)) {
            const fs::path outputDirSample = outputDirPipeline / inputSampleSingle->sampleName;

            fs::create_directories(outputDirSample);

            const fs::path outputAlignmentsPath =
                outputDirPipeline / inputSampleSingle->sampleName /
                (inputSampleSingle->sampleName + outSampleAlignedSuffix);

            samples.emplace_back(
                AlignSampleSingle{.input = *inputSampleSingle, .output = {outputAlignmentsPath}});

            const auto message = "Single-end sample " + inputSampleSingle->sampleName + " found";
            Logger::log(message);

            continue;
        }

        if (const auto* inputSamplePaired = std::get_if<AlignInputPaired>(&inputSample)) {
            const std::string parentName = inputSamplePaired->sampleName;
            const fs::path outputDirSample = outputDirPipeline / inputSamplePaired->sampleName;

            fs::create_directories(outputDirSample);

            const fs::path outputAlignmentsPath =
                outputDirSample / (parentName + outSampleAlignedSuffix);

            const fs::path outputAlignmentsMergedReadsPath =
                outputDirSample / (parentName + outSampleMergedAlignedSuffix);
            const fs::path outputAlignmentsSingletonForwardPath =
                outputDirSample / (parentName + outSampleSingletonForwardAlignedSuffix);
            const fs::path outputAlignmentsSingletonReversePath =
                outputDirSample / (parentName + outSampleSingletonReverseAlignedSuffix);
            const fs::path outputAlignmentsPairedPath =
                outputDirSample / (parentName + outSamplePairedAlignedSuffix);

            samples.emplace_back(AlignSampleMergedPaired{
                .input = *inputSamplePaired,
                .output = {.outputAlignmentsPath = outputAlignmentsPath,
                           .outputAlignmentsMergedReadsPath = outputAlignmentsMergedReadsPath,
                           .outputAlignmentsSingletonForwardReadsPath =
                               outputAlignmentsSingletonForwardPath,
                           .outputAlignmentsSingletonReverseReadsPath =
                               outputAlignmentsSingletonReversePath,
                           .outputAlignmentsPairedReadsPath = outputAlignmentsPairedPath}});

            const auto message = "Paired-end sample " + parentName + " found";
            Logger::log(message);

            continue;
        }
    }

    return samples;
}

auto AlignData::retrieveInputSamples(const fs::path& parentDir) -> std::vector<InputSampleType> {
    const std::vector<fs::path> sampleDirs = getSubDirectories(parentDir);

    std::vector<InputSampleType> samples;

    samples.reserve(sampleDirs.size());
    for (const fs::path& sampleDir : sampleDirs) {
        samples.push_back(retrieveInputSample(sampleDir));
    }

    return samples;
}

auto AlignData::retrieveInputSample(const fs::path& sampleDir) -> InputSampleType {
    const std::vector<fs::path> sampleFiles = getValidFilePaths(sampleDir, validInputSuffixes);
    const auto numSamples = sampleFiles.size();

    using namespace std::string_literals;
    if (numSamples != validInputSuffixSingleton.size() &&
        numSamples != validInputSuffixesPaired.size()) {
        const std::string message =
            "Found invalid number of sample files (" + std::to_string(numSamples) + ")" +
            ". Expectected either " + std::to_string(validInputSuffixSingleton.size()) + " or " +
            std::to_string(validInputSuffixesPaired.size()) + " in " + sampleDir.string();
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(message);
        throw std::runtime_error(message);
    }

    const auto& firstFile = sampleFiles.front();
    const auto sampleName = getSampleName(firstFile);

    // First check paired because singleton suffix is subset of paired suffix
    if (numSamples == validInputSuffixesPaired.size()) {
        return retrieveInputPaired(sampleName, sampleFiles);
    }

    return retrieveInputSingle(sampleName, firstFile);
}

auto AlignData::retrieveInputPaired(const std::string& sampleName,
                                    const std::vector<fs::path>& inputSamples) -> AlignInputPaired {
    std::unordered_map<std::string, std::optional<fs::path>> suffixToPath = {
        {preprocess::outSampleFastqPairedMergeSuffix, std::nullopt},
        {preprocess::outSampleFastqPairedForwardSingletonSuffix, std::nullopt},
        {preprocess::outSampleFastqPairedReverseSingletonSuffix, std::nullopt},
        {preprocess::outSampleFastqPairedForwardPairedSuffix, std::nullopt},
        {preprocess::outSampleFastqPairedReversePairedSuffix, std::nullopt}};

    for (const fs::path& path : inputSamples) {
        for (auto& [suffix, optPath] : suffixToPath) {
            if (hasSuffix(path, suffix)) {
                optPath = path;
                break;
            }
        }
    }

    for (const auto& [suffix, optPath] : suffixToPath) {
        if (!optPath.has_value()) {
            std::string message = "The directory " + inputSamples.front().parent_path().string();
            message += " is missing the following file: ";
            message += suffix;
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(message);
            throw std::runtime_error(message);
        }
    }

    return AlignInputPaired{
        sampleName,
        suffixToPath[preprocess::outSampleFastqPairedMergeSuffix].value(),
        suffixToPath[preprocess::outSampleFastqPairedForwardSingletonSuffix].value(),
        suffixToPath[preprocess::outSampleFastqPairedReverseSingletonSuffix].value(),
        suffixToPath[preprocess::outSampleFastqPairedForwardPairedSuffix].value(),
        suffixToPath[preprocess::outSampleFastqPairedReversePairedSuffix].value()};
}

auto AlignData::retrieveInputSingle(const std::string& sampleName, const fs::path& inputSample)
    -> AlignInputSingle {
    if (!hasSuffix(inputSample, preprocess::outSampleFastqSuffix)) {
        const std::string message = "The directory " + inputSample.parent_path().string() +
                                    " is missing the file: " + preprocess::outSampleFastqSuffix;
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(message);
        throw std::runtime_error(message);
    }
    return AlignInputSingle{sampleName, inputSample};
}

}  // namespace pipelines::align
