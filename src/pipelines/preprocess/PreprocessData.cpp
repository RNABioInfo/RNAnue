#include "PreprocessData.hpp"

#include <cstdlib>
#include <filesystem>

#include "Logger.hpp"
#include "Utility.hpp"
#include "VariantOverload.hpp"

using namespace helper;

namespace pipelines::preprocess {
auto PreprocessData::retrieveSamples(const std::string& sampleGroup, const fs::path& parentDir,
                                     const fs::path& outputDir)
    -> std::vector<PreprocessSampleType> {
    const std::vector<InputSampleType> inputSamples = retrieveInputSamples(parentDir);

    std::vector<PreprocessSampleType> samples;
    samples.reserve(inputSamples.size());

    const fs::path outputDirPipeline = outputDir / pipelinePrefix / sampleGroup;

    for (const InputSampleType& inputSample : inputSamples) {
        std::visit(
            overloaded{
                [&outputDirPipeline,
                 &samples](const PreprocessSampleInputSingle& inputSampleSingle) {
                    const std::string parentName = inputSampleSingle.sampleName;
                    const fs::path outputDirSample = outputDirPipeline / parentName;

                    fs::create_directories(outputDirSample);

                    const fs::path outputSampleFastqPath =
                        outputDirSample / (parentName + outSampleFastqSuffix);

                    const fs::path outputSampleTmpFastqDir =
                        outputDirSample / outSampleTmpFastqDirPrefix;
                    fs::create_directories(outputSampleTmpFastqDir);

                    samples.emplace_back(
                        PreprocessSampleSingle{.input = inputSampleSingle,
                                               .output = {
                                                   .tmpFastqDir = outputSampleTmpFastqDir,
                                                   .outputFastqPath = outputSampleFastqPath,
                                               }});

                    const auto message =
                        "Single-end sample " + inputSampleSingle.sampleName + " found";
                    Logger::log(message);
                },
                [&outputDirPipeline,
                 &samples](const PreprocessSampleInputPaired& inputSamplePaired) {
                    const std::string parentName = inputSamplePaired.sampleName;
                    const fs::path outputDirSample = outputDirPipeline / parentName;

                    fs::create_directories(outputDirSample);

                    const fs::path outputSampleFastqPathMerged =
                        outputDirSample / (parentName + outSampleFastqPairedMergeSuffix);
                    const fs::path outputSampleFastqPathForwardSingleton =
                        outputDirSample / (parentName + outSampleFastqPairedForwardSingletonSuffix);
                    const fs::path outputSampleFastqPathReverseSingleton =
                        outputDirSample / (parentName + outSampleFastqPairedReverseSingletonSuffix);
                    const fs::path outputSampleFastqPathForwardPaired =
                        outputDirSample / (parentName + outSampleFastqPairedForwardPairedSuffix);
                    const fs::path outputSampleFastqPathReversePaired =
                        outputDirSample / (parentName + outSampleFastqPairedReversePairedSuffix);

                    const fs::path outputSampleTmpFastqMergedDir =
                        outputDirSample / outSampleTmpMergedFastqDirPrefix;
                    fs::create_directories(outputSampleTmpFastqMergedDir);

                    const fs::path outputSampleTmpFastqSingletonForwardDir =
                        outputDirSample / outSampleTmpForwardSingletonFastqDirPrefix;
                    fs::create_directories(outputSampleTmpFastqSingletonForwardDir);

                    const fs::path outputSampleTmpFastqSingletonReverseDir =
                        outputDirSample / outSampleTmpReverseSingletonFastqDirPrefix;
                    fs::create_directories(outputSampleTmpFastqSingletonReverseDir);

                    const fs::path outputSampleTmpFastqForwardPairedDir =
                        outputDirSample / outSampleTmpForwardPairedFastqDirPrefix;
                    fs::create_directories(outputSampleTmpFastqForwardPairedDir);

                    const fs::path outputSampleTmpFastqReversePairedDir =
                        outputDirSample / outSampleTmpReversePairedFastqDirPrefix;
                    fs::create_directories(outputSampleTmpFastqReversePairedDir);

                    samples.emplace_back(PreprocessSamplePaired{
                        .input = inputSamplePaired,
                        .output = {
                            .tmpMergedFastqDir = outputSampleTmpFastqMergedDir,
                            .tmpSingletonForwardFastqDir = outputSampleTmpFastqSingletonForwardDir,
                            .tmpSingletonReverseFastqDir = outputSampleTmpFastqSingletonReverseDir,
                            .tmpPairedForwardFastqDir = outputSampleTmpFastqForwardPairedDir,
                            .tmpPairedReverseFastqDir = outputSampleTmpFastqReversePairedDir,
                            .outputMergedFastqPath = outputSampleFastqPathMerged,
                            .outputSingletonForwardFastqPath =
                                outputSampleFastqPathForwardSingleton,
                            .outputSingletonReverseFastqPath =
                                outputSampleFastqPathReverseSingleton,
                            .outputPairedForwardFastqPath = outputSampleFastqPathForwardPaired,
                            .outputPairedReverseFastqPath = outputSampleFastqPathReversePaired,
                        }});

                    const auto message = "Paired-end sample " + parentName + " found";
                    Logger::log(message);
                }},
            inputSample);
    }

    return samples;
}

auto PreprocessData::retrieveInputSamples(const fs::path& parentDir)
    -> std::vector<InputSampleType> {
    const std::vector<fs::path> sampleDirs = getSubDirectories(parentDir);

    std::vector<InputSampleType> samples;

    samples.reserve(sampleDirs.size());
    for (const fs::path& sampleDir : sampleDirs) {
        samples.push_back(retrieveInputSample(sampleDir));
    }

    return samples;
}

auto PreprocessData::retrieveInputSample(const fs::path& sampleDir) -> InputSampleType {
    const std::vector<fs::path> sampleFiles = getValidFilePaths(sampleDir, validInputSuffixes);
    const auto numSamples = sampleFiles.size();

    using namespace std::string_literals;
    if (numSamples == 0 || numSamples > 2) {
        const std::string message = (numSamples == 0 ? "No"s : "More than two"s) +
                                    " valid sample files found in " + sampleDir.string();
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(message);
        throw std::runtime_error(message);
    }

    const auto& firstFile = sampleFiles.front();
    const auto sampleName = getSampleName(firstFile);

    if (numSamples == 1) {
        if (hasAnySuffix(firstFile.stem(), validForwardTags) ||
            hasAnySuffix(firstFile.stem(), validReverseTags)) {
            const auto message =
                "Found single valid file, but expected two based on file ending in " +
                sampleDir.string();
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(message);
            throw std::runtime_error(message);
        }
        return PreprocessSampleInputSingle{sampleName, firstFile};
    }

    std::pair<fs::path, fs::path> pairedInputPaths{firstFile, sampleFiles.back()};

    if (!validatePairedFilePaths(pairedInputPaths)) {
        const auto message = "Invalid paired file paths found in " + sampleDir.string();
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(message);
        throw std::runtime_error(message);
    }

    return PreprocessSampleInputPaired{sampleName, pairedInputPaths.first, pairedInputPaths.second};
}

auto PreprocessData::validatePairedFilePaths(std::pair<fs::path, fs::path>& pairedInputPaths)
    -> bool {
    auto withoutExtensions = [](const fs::path& path) -> std::string {
        const auto stem = path.stem();
        return stem.string().substr(0, stem.string().find_last_of('.'));
    };

    const bool validForward =
        hasAnySuffix(withoutExtensions(pairedInputPaths.first), validForwardTags);

    if (!validForward) {
        pairedInputPaths = {pairedInputPaths.second, pairedInputPaths.first};

        const bool validForward =
            hasAnySuffix(withoutExtensions(pairedInputPaths.first), validForwardTags);

        if (!validForward) {
            return false;
        }
    }

    const bool validReverse =
        hasAnySuffix(withoutExtensions(pairedInputPaths.second), validReverseTags);

    return validReverse;
}
}  // namespace pipelines::preprocess
