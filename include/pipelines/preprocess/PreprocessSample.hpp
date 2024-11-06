#pragma once

// Standard
#include <filesystem>
#include <utility>
#include <variant>

namespace pipelines::preprocess {

namespace fs = std::filesystem;

struct PreprocessInput {
    std::string sampleName;
};

struct PreprocessSampleInputSingle : public PreprocessInput {
    fs::path inputFastqPath;

    PreprocessSampleInputSingle(const std::string& sampleName, fs::path samplePath)
        : PreprocessInput{sampleName}, inputFastqPath{std::move(samplePath)} {};
};

struct PreprocessSampleInputPaired : public PreprocessInput {
    fs::path inputForwardFastqPath;
    fs::path inputReverseFastqPath;

    PreprocessSampleInputPaired(const std::string& sampleName, fs::path forwardPath,
                                fs::path reversePath)
        : PreprocessInput{sampleName},
          inputForwardFastqPath{std::move(forwardPath)},
          inputReverseFastqPath{std::move(reversePath)} {};
};

using InputSampleType = std::variant<PreprocessSampleInputSingle, PreprocessSampleInputPaired>;

struct PrepocessSampleOutputSingle {
    fs::path tmpFastqDir;
    fs::path outputFastqPath;
};

struct PrepocessSampleOutputPaired {
    fs::path tmpMergedFastqDir;
    fs::path tmpSingletonForwardFastqDir;
    fs::path tmpSingletonReverseFastqDir;
    fs::path tmpPairedForwardFastqDir;
    fs::path tmpPairedReverseFastqDir;

    fs::path outputMergedFastqPath;
    fs::path outputSingletonForwardFastqPath;
    fs::path outputSingletonReverseFastqPath;
    fs::path outputPairedForwardFastqPath;
    fs::path outputPairedReverseFastqPath;
};

using OutputSampleType = std::variant<PrepocessSampleOutputSingle, PrepocessSampleOutputPaired>;

struct PreprocessSampleSingle {
    PreprocessSampleInputSingle input;
    PrepocessSampleOutputSingle output;
};

struct PreprocessSamplePaired {
    PreprocessSampleInputPaired input;
    PrepocessSampleOutputPaired output;
};

using PreprocessSampleType = std::variant<PreprocessSampleSingle, PreprocessSamplePaired>;

}  // namespace pipelines::preprocess
