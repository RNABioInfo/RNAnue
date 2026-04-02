#pragma once

// Standard
#include <string>
#include <vector>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>

// Internal
#include "GenomicFeatureGroup.hpp"

namespace dataTypes {

struct MaskedFeatureCluster {
    GenomicFeatureGroup baseFeatureGroup;
    seqan3::dna5_vector baseSequence;
    std::vector<GenomicFeatureGroup> subFeatureGroups;

    [[nodiscard]] constexpr auto isMultiCopy() const noexcept -> bool {
        return subFeatureGroups.size() > 0;
    }

    [[nodiscard]] auto csvSubFeatureRootIDs() const -> std::string {
        std::string result;
        for (const auto& group : subFeatureGroups) {
            result += group.getRoot().feature.getID() + ",";
        }
        if (!result.empty()) {
            result.pop_back();  // Remove trailing comma
        }
        return result;
    }
};

}  // namespace dataTypes
