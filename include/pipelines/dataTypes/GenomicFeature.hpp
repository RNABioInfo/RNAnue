#pragma once

// Standard
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

// Internal
#include "GenomicStrand.hpp"

namespace dataTypes {
struct GenomicFeature {
    std::string referenceID;
    std::string type;
    int32_t startPosition;
    int32_t endPosition;
    dataTypes::GenomicStrand strand;
    std::string id;
    std::optional<std::string> groupID;
    std::optional<std::string> geneName;
};

using FeatureMap = std::unordered_map<std::string, std::vector<GenomicFeature>>;

}  // namespace dataTypes
