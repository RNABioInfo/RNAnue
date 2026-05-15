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

template <typename Parameters>
inline auto getFile(const Parameters& params) -> std::filesystem::path {
    if (!params.maskMultiCopyGenes) {
        return params.referenceGenome;
    }

    auto maskedReferenceGenomeFile = maskedFilePath(params);

    if (!std::filesystem::exists(maskedReferenceGenomeFile)) {
        return params.referenceGenome;
    }

    return maskedReferenceGenomeFile;
}

}  // namespace utility::ReferenceGenomeFilePicker
