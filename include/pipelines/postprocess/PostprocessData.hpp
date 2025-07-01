#pragma once

// Standard
#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// Internal
#include "AnalyzeData.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "PipelineData.hpp"
#include "PostprocessSample.hpp"
#include "Utility.hpp"

namespace pipelines::postprocess {
namespace fs = std::filesystem;

static const std::string validInputInteractionsSuffix = analyze::outInteractionsSuffix;

static const std::array validSuffices{validInputInteractionsSuffix};

static const std::string outSuperInteractionsBEDPESuffix = "_super_interaction_regions.bedpe";
static const std::string outSuperInteractionTranscriptCountsGCSSuffix =
    "_super_interaction_transcript_counts.gct";

static const std::string pipelinePrefix = "05_postprocess";

static constexpr std::string_view conditionID = "condition";
static constexpr std::string_view controlID = "control";

struct PostprocessData : public PipelineData {
    using SampleByConditionMap = std::unordered_map<std::string, std::vector<PostprocessSample>>;

    SampleByConditionMap samples;
    fs::path superInteractionsBEDPEOutPath;
    fs::path superInteractionsCountsGCTOutPath;

    PostprocessData(const fs::path& outputDir, const fs::path& treatmentDir,
                    const std::optional<fs::path>& controlDir)
        : samples(retrieveSamples(treatmentDir, controlDir, outputDir)) {
        const fs::path outputDirPipeline = outputDir / pipelinePrefix;
        superInteractionsBEDPEOutPath =
            outputDirPipeline / ("complete" + outSuperInteractionsBEDPESuffix);
        superInteractionsCountsGCTOutPath =
            outputDirPipeline / ("complete" + outSuperInteractionTranscriptCountsGCSSuffix);
    }

   private:
    [[nodiscard]] static auto retrieveSamples(const fs::path& treatmentDir,
                                              const std::optional<fs::path>& controlDir,
                                              const fs::path& outputDir) -> SampleByConditionMap {
        const fs::path outputDirPipeline = outputDir / pipelinePrefix;

        if (fs::exists(outputDirPipeline)) {
            Logger::log<LogLevel::WARNING>(
                "Output directory already exists. Results will be overwritten.");
        }

        fs::create_directories(outputDirPipeline);

        std::unordered_map<std::string, std::vector<PostprocessSample>> samplesByCondition;

        samplesByCondition.emplace(std::string{conditionID}, retrieveSamplesInDir(treatmentDir));

        if (controlDir) {
            samplesByCondition.emplace(std::string{controlID}, retrieveSamplesInDir(*controlDir));
        }

        return samplesByCondition;
    }

    [[nodiscard]] static auto retrieveSamplesInDir(const fs::path& parentDir)
        -> std::vector<PostprocessSample> {
        const std::vector<PostprocessInput> inputSamples = retrieveInputSamples(parentDir);

        std::vector<PostprocessSample> samples;
        samples.reserve(inputSamples.size());

        for (const PostprocessInput& inputSample : inputSamples) {
            samples.emplace_back(inputSample, PostprocessOutput{});
        }

        return samples;
    }

    [[nodiscard]] static auto retrieveInputSamples(const fs::path& parentDir)
        -> std::vector<PostprocessInput> {
        const std::vector<fs::path> sampleDirs = getSubDirectories(parentDir);

        std::vector<PostprocessInput> samples;

        for (const fs::path& sampleDir : sampleDirs) {
            const std::vector<fs::path> sampleFiles =
                helper::getValidFilePaths(sampleDir, validSuffices);

            if (sampleFiles.size() != validSuffices.size()) {
                const std::string message = "Expected " + std::to_string(validSuffices.size()) +
                                            " files in " + sampleDir.string() + " but found " +
                                            std::to_string(sampleFiles.size()) + " valid files";
                Logger::log<IncludeSourceLocation, LogLevel::ERROR>(message);
            }

            const std::string sampleName = getSampleName(sampleFiles.front());

            std::optional<fs::path> interactionsPath;

            for (const fs::path& sampleFile : sampleFiles) {
                if (helper::hasSuffix(sampleFile, validInputInteractionsSuffix)) {
                    interactionsPath = sampleFile;
                } else {
                    const std::string message = "Unexpected file " + sampleFile.string() +
                                                " found in " + sampleDir.string();
                    Logger::log<LogLevel::WARNING>(message);
                }
            }

            if (!interactionsPath.has_value()) {
                const std::string message =
                    "Missing one or more required files in " + sampleDir.string();
                Logger::log<IncludeSourceLocation, LogLevel::ERROR>(message);
            }

            samples.emplace_back(sampleName, interactionsPath.value());
        }

        return samples;
    }
};

}  // namespace pipelines::postprocess
