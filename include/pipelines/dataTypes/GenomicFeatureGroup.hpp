#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "GenomicFeature.hpp"

namespace dataTypes {

class GenomicFeatureGroup final {
   public:
    using NodeIndex = std::uint32_t;

    static constexpr NodeIndex invalidIndex = static_cast<NodeIndex>(-1);

    struct Node final {
        GenomicFeature feature;
        NodeIndex parentIndex = invalidIndex;

        std::uint32_t firstChildOffset = 0;
        std::uint32_t childCount = 0;
    };

    struct BuildIssue final {
        enum class Kind : std::uint8_t { duplicateId, missingParent, cycleDetected };

        Kind kind{};
        std::string id;
        std::optional<std::string> parentId;
    };

    struct BuildResult;

    [[nodiscard]] static auto buildFromFlat(std::vector<GenomicFeature>&& features,
                                            std::string_view expectedRootId = {}) -> BuildResult;

    [[nodiscard]] auto getRootIndex() const noexcept -> NodeIndex { return rootIndex; }
    [[nodiscard]] auto getRoot() const noexcept -> const Node& { return nodes.at(rootIndex); }

    [[nodiscard]] auto getNode(NodeIndex nodeIndex) const noexcept -> const Node& {
        return nodes.at(nodeIndex);
    }

    [[nodiscard]] auto children(NodeIndex nodeIndex) const noexcept -> std::span<const NodeIndex> {
        const Node& node = nodes.at(nodeIndex);
        return std::span<const NodeIndex>{childIndices}.subspan(node.firstChildOffset,
                                                                node.childCount);
    }

    [[nodiscard]] auto allChildrenOfRoot() const noexcept -> std::span<const NodeIndex> {
        const Node& node = nodes.at(rootIndex);
        return std::span<const NodeIndex>{childIndices}.subspan(node.firstChildOffset,
                                                                node.childCount);
    }

    [[nodiscard]] auto tryFindById(std::string_view identifier) const noexcept
        -> std::optional<NodeIndex> {
        const auto iterator = idToIndex.find(std::string{identifier});
        if (iterator == idToIndex.end()) {
            return std::nullopt;
        }
        return iterator->second;
    }

    template <typename Function>
    auto dfs(NodeIndex startIndex, Function&& visit) const -> void {
        std::vector<NodeIndex> stack;
        constexpr std::size_t capacity = 32;
        stack.reserve(capacity);
        stack.push_back(startIndex);

        while (!stack.empty()) {
            const NodeIndex nodeIndex = stack.back();
            stack.pop_back();

            visit(nodes.at(nodeIndex));

            const auto childSpan = children(nodeIndex);
            for (auto reverseIndex = childSpan.size(); reverseIndex > 0; --reverseIndex) {
                stack.push_back(childSpan[reverseIndex - 1]);
            }
        }
    }

    [[nodiscard]] auto getNodes() const noexcept -> const std::vector<Node>& { return nodes; }

   private:
    std::vector<Node> nodes;
    std::vector<NodeIndex> childIndices;
    std::unordered_map<std::string, NodeIndex> idToIndex;
    NodeIndex rootIndex = invalidIndex;
};

struct GenomicFeatureGroup::BuildResult final {
    GenomicFeatureGroup group;
    std::vector<BuildIssue> issues;
};

using ParentIDToFeatureGroupMap = std::unordered_map<std::string, GenomicFeatureGroup>;

}  // namespace dataTypes
