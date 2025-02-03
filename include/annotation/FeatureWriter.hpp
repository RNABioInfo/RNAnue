#pragma once

// Standard
#include <string>

// Internal
#include "FeatureAnnotator.hpp"
#include "FileType.hpp"

namespace annotation {

class FeatureWriter {
   public:
    FeatureWriter() = delete;
    FeatureWriter(const FeatureWriter &) = delete;
    FeatureWriter(FeatureWriter &&) = delete;
    auto operator=(const FeatureWriter &) -> FeatureWriter & = delete;
    auto operator=(FeatureWriter &&) -> FeatureWriter & = delete;
    ~FeatureWriter() = delete;

    static void write(const FeatureTreeMap &featureTreeMap, const std::string &outputPath,
                      FileType::Value fileType);
};

}  // namespace annotation
