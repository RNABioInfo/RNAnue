#pragma once

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_set>
#include <vector>

#include "FileType.hpp"
#include "GenomicFeature.hpp"

namespace annotation {

enum class ParsedParentOrigin : unsigned char {
    None,
    GffParent,
    GtfTranscriptId,
};

struct ParsedFeatureRecord final {
    dataTypes::GenomicFeature feature;
    std::vector<std::string> explicitParentIds;
    std::size_t lineNumber{};
    ParsedParentOrigin parentOrigin{ParsedParentOrigin::None};
};

struct HierarchyResolutionStats final {
    std::size_t inferredParents{};
    std::size_t aliasRepairs{};
    std::size_t equivalentDuplicates{};
    std::vector<std::string> repairExamples;
};

struct HierarchyResolutionResult final {
    std::vector<dataTypes::GenomicFeature> features;
    HierarchyResolutionStats stats;
};

class AnnotationHierarchyResolver final {
   public:
    [[nodiscard]] static auto resolve(std::vector<ParsedFeatureRecord>&& records,
                                      const std::filesystem::path& featureFilePath,
                                      FileType fileType,
                                      const std::unordered_set<std::string>& includedFeatures)
        -> HierarchyResolutionResult;
};

}  // namespace annotation
