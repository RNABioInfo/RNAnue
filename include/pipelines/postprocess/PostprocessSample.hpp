#pragma once

// Standard
#include <filesystem>
#include <string>

namespace pipelines::postprocess {
namespace fs = std::filesystem;

struct PostprocessInput {
    std::string sampleName;

    fs::path interactionsPath;
};

struct PostprocessOutput {
    // Currently no sample specific outputs is planned all outputs should be summarized
};

struct PostprocessSample {
    PostprocessInput input;
    PostprocessOutput output;
};

}  // namespace pipelines::postprocess
