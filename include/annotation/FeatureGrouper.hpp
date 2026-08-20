#pragma once

// Standard
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicFeatureGroup.hpp"

namespace annotation::FeatureGrouper {

using GroupingFeature = dataTypes::GenomicFeature;

using GroupMap = std::unordered_map<std::string, dataTypes::GenomicFeatureGroup>;

struct GroupingResult {
    GroupMap groups;
    std::unordered_set<std::string> groupKeys;
};

[[nodiscard]] auto groupByHierarchy(std::vector<GroupingFeature>&& groupingFeatures)
    -> GroupingResult;

[[nodiscard]] auto groupByValidatedHierarchy(std::vector<GroupingFeature>&& groupingFeatures)
    -> GroupingResult;

[[nodiscard]] auto groupByDirectParentID(std::vector<GroupingFeature>&& groupingFeatures)
    -> GroupingResult;

}  // namespace annotation::FeatureGrouper
