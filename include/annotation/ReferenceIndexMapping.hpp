// ReferenceIndexMapping.hpp
#pragma once

// Standard
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>

// Internal
#include "LogLevel.hpp"
#include "Logger.hpp"

namespace annotation {

using ReferenceIDToIndexMap = std::unordered_map<std::string, int>;

enum class MissingReferencePolicy : std::uint8_t {
    Reject,
    Create,
};

class ReferenceIndexMapping {
   public:
    explicit ReferenceIndexMapping(ReferenceIDToIndexMap initialMap,
                                   MissingReferencePolicy policyParam) noexcept
        : mapping(std::move(initialMap)),
          policy(policyParam),
          nextIndex(computeNextIndex(mapping)) {}

    [[nodiscard]] static auto defaultCreate() noexcept {
        return ReferenceIndexMapping{{}, MissingReferencePolicy::Create};
    }

    [[nodiscard]] auto findIndex(std::string_view referenceId) const noexcept
        -> std::optional<int> {
        lookupKey.assign(referenceId.data(), referenceId.size());
        const auto iter = mapping.find(lookupKey);
        if (iter == mapping.end()) {
            return std::nullopt;
        }
        return iter->second;
    }

    [[nodiscard]] auto getIndex(std::string_view referenceId) -> int {
        if (const auto found = findIndex(referenceId)) {
            return found.value();
        }

        if (policy == MissingReferencePolicy::Create) {
            return createIndex(referenceId);
        }

        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Unknown reference id: {}",
                                                            std::string(referenceId));

        std::unreachable();
    }

    [[nodiscard]] auto findID(int index) const -> std::optional<std::string_view> {
        for (const auto& [identifier, idx] : mapping) {
            if (idx == index) {
                return identifier;
            }
        }
        return std::nullopt;
    }

    [[nodiscard]] auto sortedReferenceIDs() const -> std::deque<std::string> {
        std::deque<std::pair<int, std::string_view>> indexedIDs;

        for (const auto& [identifier, idx] : mapping) {
            indexedIDs.emplace_back(idx, identifier);
        }

        std::ranges::sort(indexedIDs, {}, &std::pair<int, std::string_view>::first);

        std::deque<std::string> sortedIDs;
        sortedIDs.resize(indexedIDs.size());
        for (std::size_t i = 0; i < indexedIDs.size(); ++i) {
            sortedIDs[i] = std::string{indexedIDs[i].second};
        }

        return sortedIDs;
    }

    [[nodiscard]] auto size() const noexcept -> std::size_t { return mapping.size(); }
    [[nodiscard]] auto data() const noexcept -> const ReferenceIDToIndexMap& { return mapping; }

    auto reserve(std::size_t count) -> void { mapping.reserve(count); }

   private:
    ReferenceIDToIndexMap mapping;
    MissingReferencePolicy policy{MissingReferencePolicy::Reject};
    int nextIndex{0};

    mutable std::string lookupKey;

    static auto computeNextIndex(const ReferenceIDToIndexMap& mapData) noexcept -> int {
        if (mapData.empty()) {
            return 0;
        }

        const auto iterMax = std::ranges::max_element(
            mapData, [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });

        const long long candidate = static_cast<long long>(iterMax->second) + 1LL;
        if (candidate >= static_cast<long long>(std::numeric_limits<int>::max())) {
            return std::numeric_limits<int>::max();
        }
        return static_cast<int>(candidate);
    }

    auto createIndex(std::string_view referenceId) -> int {
        const int assigned = nextIndex;

        if (nextIndex < std::numeric_limits<int>::max()) {
            ++nextIndex;
        }

        mapping.try_emplace(std::string{referenceId}, assigned);
        return assigned;
    }
};

}  // namespace annotation
