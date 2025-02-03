#pragma once

// Standard
#include <filesystem>
#include <optional>
#include <string>

// Internal
#include "AlignParameters.hpp"
#include "AnalyzeParameters.hpp"
#include "CompleteParameters.hpp"
#include "DetectParameters.hpp"
#include "PipelineData.hpp"
#include "PreprocessParameters.hpp"
#include "Utility.hpp"

using namespace pipelines;

class Runner {
   public:
    Runner() = delete;
    Runner(const Runner&) = delete;
    Runner(Runner&&) = delete;
    auto operator=(const Runner&) -> Runner& = delete;
    auto operator=(Runner&&) -> Runner& = delete;
    ~Runner() = delete;

    static void runPipeline(int argc, const char* const argv[]);

   private:
    struct InputDirectories {
        InputDirectories(const fs::path& parentDir, const std::string& pipelinePrefix)
            : treatmentInputDir(parentDir / pipelinePrefix / treatmentSampleGroup),
              controlInputDir(
                  helper::getDirIfExists(parentDir / pipelinePrefix / controlSampleGroup)) {};

        fs::path treatmentInputDir;
        std::optional<fs::path> controlInputDir;
    };

    struct Pipeline {
        void operator()(const preprocess::PreprocessParameters& params);
        void operator()(const align::AlignParameters& params);
        void operator()(const detect::DetectParameters& params);
        void operator()(const analyze::AnalyzeParameters& params);
        void operator()(const CompleteParameters& params);
    };

    static void runPreprocessPipeline(const preprocess::PreprocessParameters& parameters);
    static void runAlignPipeline(const align::AlignParameters& parameters);
    static void runDetectPipeline(const detect::DetectParameters& parameters);
    static void runAnalyzePipeline(const analyze::AnalyzeParameters& parameters);

    static void runCompletePipeline(const CompleteParameters& parameters);
};
