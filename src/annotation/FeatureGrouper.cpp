#include "FeatureGrouper.hpp"

#include <algorithm>
#include <cstddef>
#include <format>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "GenomicFeature.hpp"
#include "GenomicFeatureGroup.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"

namespace annotation::FeatureGrouper {

namespace {

using IdentifierToIndexMap = std::unordered_map<std::string, std::size_t>;
using GroupKeyToIndicesMap = std::unordered_map<std::string, std::vector<std::size_t>>;

struct GroupingWarnings {
    std::unordered_set<std::string> warnedDuplicateIdentifiers;
    std::unordered_set<std::string> warnedCycleGroupKeys;
    std::unordered_set<std::string> warnedOrphanGroupKeys;
};

auto warnDuplicateIfNeeded(GroupingWarnings& warnings, const std::string& identifier) -> void {
    if (!warnings.warnedDuplicateIdentifiers.insert(identifier).second) {
        return;
    }
    Logger::log<LogLevel::WARNING>(std::format(
        "Duplicate feature ID encountered: '{}'. Parent resolution will use the first occurrence.",
        identifier));
}

auto warnOrphanIfNeeded(GroupingWarnings& warnings, const std::string& groupKey) -> void {
    if (!warnings.warnedOrphanGroupKeys.insert(groupKey).second) {
        return;
    }
    Logger::log<LogLevel::WARNING>(std::format(
        "Found features whose Parent '{}' is missing from the file. Grouping them under '{}'.",
        groupKey, groupKey));
}

auto warnCycleIfNeeded(GroupingWarnings& warnings, const std::string& groupKey) -> void {
    if (!warnings.warnedCycleGroupKeys.insert(groupKey).second) {
        return;
    }
    Logger::log<LogLevel::WARNING>(std::format(
        "Circular Parent relationship detected involving '{}'. Grouping cycle under '{}'.",
        groupKey, groupKey));
}

auto buildIdentifierToIndex(const std::vector<GroupingFeature>& groupingFeatures,
                            GroupingWarnings& warnings) -> IdentifierToIndexMap {
    IdentifierToIndexMap identifierToIndex;
    identifierToIndex.reserve(groupingFeatures.size());

    for (std::size_t featureIndex = 0; featureIndex < groupingFeatures.size(); ++featureIndex) {
        const std::string& identifier = groupingFeatures[featureIndex].getID();
        const auto [iterator, inserted] = identifierToIndex.try_emplace(identifier, featureIndex);
        if (!inserted) {
            warnDuplicateIfNeeded(warnings, identifier);
        }
    }

    return identifierToIndex;
}

auto resolveGroupKeyForIndex(std::size_t startIndex,
                             const std::vector<GroupingFeature>& groupingFeatures,
                             const IdentifierToIndexMap& identifierToIndex,
                             std::vector<std::string>& resolvedGroupKey, GroupingWarnings& warnings)
    -> const std::string& {
    auto& cachedKey = resolvedGroupKey[startIndex];
    if (!cachedKey.empty()) {
        return cachedKey;
    }

    std::vector<std::size_t> traversalIndices;
    constexpr std::size_t capacity = 8;
    traversalIndices.reserve(capacity);

    std::size_t currentIndex = startIndex;

    while (true) {
        auto& currentCachedKey = resolvedGroupKey[currentIndex];
        if (!currentCachedKey.empty()) {
            cachedKey = currentCachedKey;
            break;
        }

        traversalIndices.push_back(currentIndex);

        const auto& parentIdOpt = groupingFeatures[currentIndex].getParentID();
        if (!parentIdOpt) {
            cachedKey = groupingFeatures[currentIndex].getID();  // true root
            break;
        }

        const std::string& parentIdentifier = *parentIdOpt;
        const auto parentIterator = identifierToIndex.find(parentIdentifier);

        if (parentIterator == identifierToIndex.end()) {
            cachedKey = parentIdentifier;  // orphan chain
            warnOrphanIfNeeded(warnings, cachedKey);
            break;
        }

        const std::size_t parentIndex = parentIterator->second;
        const bool cycleDetected =
            std::ranges::find(traversalIndices, parentIndex) != traversalIndices.end();

        if (cycleDetected) {
            cachedKey = parentIdentifier;  // stable cycle key
            warnCycleIfNeeded(warnings, cachedKey);
            break;
        }

        currentIndex = parentIndex;
    }

    for (const std::size_t traversalIndex : traversalIndices) {
        resolvedGroupKey[traversalIndex] = cachedKey;  // path compression
    }

    return cachedKey;
}

auto bucketByResolvedGroupKey(const std::vector<GroupingFeature>& groupingFeatures,
                              const IdentifierToIndexMap& identifierToIndex,
                              GroupingWarnings& warnings) -> GroupKeyToIndicesMap {
    std::vector<std::string> resolvedGroupKey;
    resolvedGroupKey.resize(groupingFeatures.size());

    GroupKeyToIndicesMap groupKeyToIndices;
    groupKeyToIndices.reserve(groupingFeatures.size());

    for (std::size_t featureIndex = 0; featureIndex < groupingFeatures.size(); ++featureIndex) {
        const std::string& groupKey = resolveGroupKeyForIndex(
            featureIndex, groupingFeatures, identifierToIndex, resolvedGroupKey, warnings);

        groupKeyToIndices[groupKey].push_back(featureIndex);
    }

    return groupKeyToIndices;
}

auto logBuildIssues(const std::string& groupKey,
                    const std::vector<dataTypes::GenomicFeatureGroup::BuildIssue>& issues) -> void {
    using Kind = dataTypes::GenomicFeatureGroup::BuildIssue::Kind;

    for (const auto& issue : issues) {
        switch (issue.kind) {
            case Kind::duplicateId:
                Logger::log<LogLevel::WARNING>(std::format(
                    "Group '{}': duplicate feature ID '{}' encountered while building tree.",
                    groupKey, issue.id));
                break;
            case Kind::missingParent:
                Logger::log<LogLevel::WARNING>(
                    std::format("Group '{}': feature '{}' has missing Parent '{}'.", groupKey,
                                issue.id, issue.parentId.value_or(std::string{"<null>"})));
                break;
            case Kind::cycleDetected:
                Logger::log<LogLevel::WARNING>(std::format(
                    "Group '{}': cycle detected involving feature '{}' (Parent '{}').", groupKey,
                    issue.id, issue.parentId.value_or(std::string{"<null>"})));
                break;
        }
    }
}

auto buildTreeGroup(std::string_view groupKey, const std::vector<std::size_t>& featureIndices,
                    std::vector<GroupingFeature>& groupingFeatures)
    -> dataTypes::GenomicFeatureGroup {
    std::vector<dataTypes::GenomicFeature> groupFeatures;
    groupFeatures.reserve(featureIndices.size());

    for (const std::size_t featureIndex : featureIndices) {
        groupFeatures.emplace_back(groupingFeatures[featureIndex]);
    }

    auto buildResult =
        dataTypes::GenomicFeatureGroup::buildFromFlat(std::move(groupFeatures), groupKey);
    logBuildIssues(std::string{groupKey}, buildResult.issues);

    return std::move(buildResult.group);
}

auto copyWithoutParent(const GroupingFeature& feature) -> dataTypes::GenomicFeature {
    return dataTypes::GenomicFeature{feature.getType(),
                                     feature.getGenomicRegion(),
                                     feature.getID(),
                                     std::nullopt,
                                     feature.getGeneName(),
                                     feature.getAttributes()};
}

auto buildDirectParentGroup(std::string_view groupKey,
                            const std::vector<std::size_t>& featureIndices,
                            const std::vector<GroupingFeature>& groupingFeatures)
    -> dataTypes::GenomicFeatureGroup {
    std::vector<dataTypes::GenomicFeature> groupFeatures;
    groupFeatures.reserve(featureIndices.size());

    const bool groupContainsParent =
        std::ranges::any_of(featureIndices, [&](const std::size_t featureIndex) {
            return groupingFeatures[featureIndex].getID() == groupKey;
        });

    for (const std::size_t featureIndex : featureIndices) {
        const auto& feature = groupingFeatures[featureIndex];
        if (groupContainsParent) {
            groupFeatures.emplace_back(feature);
        } else {
            groupFeatures.emplace_back(copyWithoutParent(feature));
        }
    }

    auto buildResult =
        dataTypes::GenomicFeatureGroup::buildFromFlat(std::move(groupFeatures), groupKey);
    logBuildIssues(std::string{groupKey}, buildResult.issues);

    return std::move(buildResult.group);
}

auto directParentGroupKey(const GroupingFeature& feature) -> const std::string& {
    const auto& parentId = feature.getParentID();
    if (parentId && !parentId->empty()) {
        return *parentId;
    }

    return feature.getID();
}

}  // namespace

auto groupByHierarchy(std::vector<GroupingFeature>&& groupingFeatures) -> GroupingResult {
    GroupingWarnings warnings{};

    const IdentifierToIndexMap identifierToIndex =
        buildIdentifierToIndex(groupingFeatures, warnings);

    GroupKeyToIndicesMap groupKeyToIndices =
        bucketByResolvedGroupKey(groupingFeatures, identifierToIndex, warnings);

    GroupingResult result{};
    result.groups.reserve(groupKeyToIndices.size());
    result.groupKeys.reserve(groupKeyToIndices.size());

    for (auto& [groupKey, featureIndices] : groupKeyToIndices) {
        result.groupKeys.insert(groupKey);
        result.groups.try_emplace(groupKey,
                                  buildTreeGroup(groupKey, featureIndices, groupingFeatures));
    }

    return result;
}

auto groupByDirectParentID(std::vector<GroupingFeature>&& groupingFeatures) -> GroupingResult {
    GroupKeyToIndicesMap groupKeyToIndices;
    groupKeyToIndices.reserve(groupingFeatures.size());

    for (std::size_t featureIndex = 0; featureIndex < groupingFeatures.size(); ++featureIndex) {
        groupKeyToIndices[directParentGroupKey(groupingFeatures[featureIndex])].push_back(
            featureIndex);
    }

    GroupingResult result{};
    result.groups.reserve(groupKeyToIndices.size());
    result.groupKeys.reserve(groupKeyToIndices.size());

    for (auto& [groupKey, featureIndices] : groupKeyToIndices) {
        result.groupKeys.insert(groupKey);
        result.groups.try_emplace(
            groupKey, buildDirectParentGroup(groupKey, featureIndices, groupingFeatures));
    }

    return result;
}

}  // namespace annotation::FeatureGrouper
