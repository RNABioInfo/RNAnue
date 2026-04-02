#pragma once

#include <filesystem>

#include "AlignData.hpp"
#include "Constants.hpp"
#include "GeneralParameters.hpp"

namespace utility::ReferenceGenomeFilePicker {

inline auto maskedFilePath(const GeneralParameters& params) -> std::filesystem::path {
    auto alignOutputDir = params.outputDir / pipelines::align::pipelinePrefix;
    return alignOutputDir / constants::pipelines::maskedReferenceGenomeFileName;
}

inline auto getFile(const GeneralParameters& params) -> std::filesystem::path {
    // Check if masked annotation file exists in align output
    auto alignOutputDir = params.outputDir / pipelines::align::pipelinePrefix;

    auto maskedReferenceGenomeFile =
        alignOutputDir / constants::pipelines::maskedReferenceGenomeFileName;

    if (!std::filesystem::exists(maskedReferenceGenomeFile)) {
        return params.featuresInPath;
    }

    return maskedReferenceGenomeFile;
}

}  // namespace utility::ReferenceGenomeFilePicker
