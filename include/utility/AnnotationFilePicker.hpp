#pragma once

#include <filesystem>

#include "AlignData.hpp"
#include "Constants.hpp"
#include "GeneralParameters.hpp"

namespace utility::AnnotationFilePicker {

inline auto maskedFilePath(const GeneralParameters& params) -> std::filesystem::path {
    auto alignOutputDir = params.outputDir / pipelines::align::pipelinePrefix;
    return alignOutputDir / constants::pipelines::maskedAnnotationFileName;
}

inline auto getFile(const GeneralParameters& params) -> std::filesystem::path {
    // Check if masked annotation file exists in align output
    auto maskedAnnotationFile = maskedFilePath(params);

    if (!std::filesystem::exists(maskedAnnotationFile)) {
        return params.featuresInPath;
    }

    return maskedAnnotationFile;
}

}  // namespace utility::AnnotationFilePicker
