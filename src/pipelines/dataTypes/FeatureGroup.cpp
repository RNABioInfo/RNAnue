#include <cstdint>
#include <format>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "GenomicFeature.hpp"
#include "GenomicFeatureGroup.hpp"

namespace dataTypes {

namespace {

auto recordIssue(std::vector<GenomicFeatureGroup::BuildIssue>& issues,
                 GenomicFeatureGroup::BuildIssue::Kind kind, std::string identifier,
                 std::optional<std::string> parentId) -> void {
    issues.emplace_back(GenomicFeatureGroup::BuildIssue{
        .kind = kind,
        .id = std::move(identifier),
        .parentId = std::move(parentId),
    });
}

auto chooseRootIndex(
    const std::vector<GenomicFeatureGroup::Node>& nodes, std::string_view expectedRootId,
    const std::unordered_map<std::string, GenomicFeatureGroup::NodeIndex>& idToIndex)
    -> GenomicFeatureGroup::NodeIndex {
    if (!expectedRootId.empty()) {
        const auto iterator = idToIndex.find(std::string{expectedRootId});
        if (iterator != idToIndex.end()) {
            return iterator->second;
        }
    }

    // Prefer a node without parent (typical GFF root).
    for (GenomicFeatureGroup::NodeIndex nodeIndex = 0; nodeIndex < nodes.size(); ++nodeIndex) {
        if (nodes[nodeIndex].parentIndex == GenomicFeatureGroup::invalidIndex) {
            return nodeIndex;
        }
    }

    // Fallback to first node.
    return nodes.empty() ? GenomicFeatureGroup::invalidIndex : 0;
}

}  // namespace

auto GenomicFeatureGroup::requireRootedTree(std::string_view context) const -> void {
    auto fail = [&](std::string_view detail) {
        throw std::invalid_argument{std::format("{}: invalid feature hierarchy: {}", context,
                                                detail)};
    };
    if (rootIndex >= nodes.size()) {
        fail("missing root");
    }
    if (nodes[rootIndex].parentIndex != invalidIndex ||
        nodes[rootIndex].feature.getParentID().has_value()) {
        fail(std::format("root '{}' still declares a Parent",
                         nodes[rootIndex].feature.getID()));
    }
    for (NodeIndex index = 0; index < nodes.size(); ++index) {
        if (index == rootIndex) {
            continue;
        }
        const auto& node = nodes[index];
        if (node.parentIndex >= nodes.size()) {
            fail(std::format("feature '{}' has no resolved Parent (requested '{}')",
                             node.feature.getID(), node.feature.getParentID().value_or("<none>")));
        }
    }

    // Follow each parent chain once, including disconnected cycles. No recursive traversal.
    std::vector<std::uint8_t> state(nodes.size(), 0);
    state[rootIndex] = 2;
    std::vector<NodeIndex> path;
    for (NodeIndex index = 0; index < nodes.size(); ++index) {
        path.clear();
        auto current = index;
        while (state[current] == 0) {
            state[current] = 1;
            path.push_back(current);
            current = nodes[current].parentIndex;
        }
        if (state[current] == 1) {
            fail(std::format("cycle involving feature '{}'", nodes[current].feature.getID()));
        }
        for (const auto visited : path) {
            state[visited] = 2;
        }
    }
}

auto GenomicFeatureGroup::buildFromFlat(std::vector<GenomicFeature>&& features,
                                        std::string_view expectedRootId) -> BuildResult {
    BuildResult result{};
    GenomicFeatureGroup& group = result.group;

    if (features.empty()) {
        return result;
    }

    group.nodes.reserve(features.size());
    group.idToIndex.reserve(features.size());

    // Create nodes and ID index (first wins).
    for (auto& feature : features) {
        const std::string& identifier = feature.getID();

        const auto [iterator, inserted] =
            group.idToIndex.try_emplace(identifier, static_cast<NodeIndex>(group.nodes.size()));

        if (!inserted) {
            recordIssue(result.issues, BuildIssue::Kind::duplicateId, identifier,
                        feature.getParentID());
            continue;
        }

        group.nodes.emplace_back(Node{
            .feature = std::move(feature),
            .parentIndex = invalidIndex,
            .firstChildOffset = 0,
            .childCount = 0,
        });
    }

    // Resolve parent indices (with basic cycle detection).
    // We treat parentId == feature.getGroupID() (your semantics) as "parent one level up".
    std::vector<std::uint8_t> visitState(group.nodes.size(), 0);  // 0=unvisited,1=visiting,2=done

    auto resolveParent = [&](NodeIndex nodeIndex, auto&& resolveParentRef) -> void {
        if (visitState[nodeIndex] == 2) {
            return;
        }
        if (visitState[nodeIndex] == 1) {
            // Cycle detected: break by clearing parent.
            recordIssue(result.issues, BuildIssue::Kind::cycleDetected,
                        group.nodes[nodeIndex].feature.getID(),
                        group.nodes[nodeIndex].feature.getParentID());
            group.nodes[nodeIndex].parentIndex = invalidIndex;
            visitState[nodeIndex] = 2;
            return;
        }

        visitState[nodeIndex] = 1;

        const auto& parentIdOpt = group.nodes[nodeIndex].feature.getParentID();
        if (!parentIdOpt) {
            group.nodes[nodeIndex].parentIndex = invalidIndex;
            visitState[nodeIndex] = 2;
            return;
        }

        const std::string& parentId = *parentIdOpt;
        const auto parentIterator = group.idToIndex.find(parentId);

        if (parentIterator == group.idToIndex.end()) {
            recordIssue(result.issues, BuildIssue::Kind::missingParent,
                        group.nodes[nodeIndex].feature.getID(), parentId);
            group.nodes[nodeIndex].parentIndex = invalidIndex;
            visitState[nodeIndex] = 2;
            return;
        }

        const NodeIndex parentIndex = parentIterator->second;
        resolveParentRef(parentIndex, resolveParentRef);
        group.nodes[nodeIndex].parentIndex = parentIndex;

        visitState[nodeIndex] = 2;
    };

    for (NodeIndex nodeIndex = 0; nodeIndex < group.nodes.size(); ++nodeIndex) {
        resolveParent(nodeIndex, resolveParent);
    }

    group.rootIndex = chooseRootIndex(group.nodes, expectedRootId, group.idToIndex);
    if (group.rootIndex == invalidIndex) {
        return result;
    }

    // Count children.
    std::vector<std::uint32_t> childCounts(group.nodes.size(), 0);
    for (NodeIndex nodeIndex = 0; nodeIndex < group.nodes.size(); ++nodeIndex) {
        const NodeIndex parentIndex = group.nodes[nodeIndex].parentIndex;
        if (parentIndex != invalidIndex) {
            ++childCounts[parentIndex];
        }
    }

    // Prefix sums to compute offsets.
    std::uint32_t runningOffset = 0;
    for (NodeIndex nodeIndex = 0; nodeIndex < group.nodes.size(); ++nodeIndex) {
        group.nodes[nodeIndex].firstChildOffset = runningOffset;
        group.nodes[nodeIndex].childCount = childCounts[nodeIndex];
        runningOffset += childCounts[nodeIndex];
    }

    group.childIndices.assign(runningOffset, invalidIndex);

    // Fill child indices using a write cursor per parent.
    std::vector<std::uint32_t> writeCursor(group.nodes.size(), 0);
    for (NodeIndex nodeIndex = 0; nodeIndex < group.nodes.size(); ++nodeIndex) {
        const NodeIndex parentIndex = group.nodes[nodeIndex].parentIndex;
        if (parentIndex == invalidIndex) {
            continue;
        }

        const std::uint32_t offset =
            group.nodes[parentIndex].firstChildOffset + writeCursor[parentIndex]++;
        group.childIndices[offset] = nodeIndex;
    }

    return result;
}

}  // namespace dataTypes
