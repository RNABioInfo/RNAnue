#pragma once

// Standard
#include <algorithm>
#include <cassert>
#include <string>
#include <utility>
#include <vector>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "Region.hpp"

namespace dataTypes {

/**
 * @brief Represents a group of genomic features that share the same group ID, strand, and reference
 * ID index.
 *
 * The genomic region of the group is the merged region of all features.
 */
struct GenomicFeatureGroup {
    /**
     * @brief Constructs a GenomicFeatureGroup from a vector of GenomicFeature objects.
     *
     * @param features A vector of features; must not be empty and all features must share the same
     * group ID, strand, and reference ID index.
     */
    explicit GenomicFeatureGroup(std::vector<GenomicFeature> features)
        : groupID{getGroupID(features)},
          genomicRegion{getReferenceIDIndex(features), getMergedRegion(features),
                        getStrand(features)},
          features{std::move(features)} {}

    std::string groupID;                   ///< Group identifier shared by all features.
    GenomicRegion genomicRegion;           ///< Merged genomic region of all features.
    std::vector<GenomicFeature> features;  ///< The list of individual genomic features.

   private:
    /**
     * @brief Retrieves and validates the common group ID from the features.
     *
     * @param features A vector of genomic features.
     * @return The common group ID.
     */
    [[nodiscard]] static auto getGroupID(const std::vector<GenomicFeature>& features) noexcept
        -> std::string {
        assert(!features.empty() && "Expected at least one feature to create a feature group.");
        const auto& firstFeature = features.front();
        assert(firstFeature.getGroupID().has_value() &&
               "All features must have a groupID to be coalesced into a GenomicFeatureGroup.");
        const std::string groupID = firstFeature.getGroupID().value();

        const bool allSame =
            std::ranges::all_of(features, [&groupID](const GenomicFeature& feature) {
                return feature.getGroupID().value_or("") == groupID;
            });
        assert(
            allSame &&
            "All features must have the same groupID to be coalesced into a GenomicFeatureGroup.");
        return groupID;
    }

    /**
     * @brief Merges the genomic regions of all features.
     *
     * @param features A vector of genomic features.
     * @return The merged region.
     */
    [[nodiscard]] static auto getMergedRegion(const std::vector<GenomicFeature>& features) noexcept
        -> Region {
        assert(!features.empty() && "Expected at least one feature to create a feature group.");
        Region mergedRegion = features.front().getGenomicRegion().getRegion();
        for (const auto& feature : features) {
            mergedRegion.merge(feature.getGenomicRegion().getRegion());
        }
        return mergedRegion;
    }

    /**
     * @brief Retrieves and validates the common strand from the features.
     *
     * @param features A vector of genomic features.
     * @return The common GenomicStrand.
     */
    [[nodiscard]] static auto getStrand(const std::vector<GenomicFeature>& features) noexcept
        -> GenomicStrand {
        assert(!features.empty() && "Expected at least one feature to create a feature group.");
        const auto strand = features.front().getGenomicRegion().getStrand();
        const bool allSame = std::ranges::all_of(features, [strand](const GenomicFeature& feature) {
            return feature.getGenomicRegion().getStrand() == strand;
        });
        assert(
            allSame &&
            "All features must have the same strand to be coalesced into a GenomicFeatureGroup.");
        return strand;
    }

    /**
     * @brief Retrieves and validates the common reference ID index from the features.
     *
     * @param features A vector of genomic features.
     * @return The common reference ID index.
     */
    [[nodiscard]] static auto getReferenceIDIndex(
        const std::vector<GenomicFeature>& features) noexcept -> int {
        assert(!features.empty() && "Expected at least one feature to create a feature group.");
        const int refIDIndex = features.front().getGenomicRegion().getReferenceIDIndex();
        const bool allSame =
            std::ranges::all_of(features, [refIDIndex](const GenomicFeature& feature) {
                return feature.getGenomicRegion().getReferenceIDIndex() == refIDIndex;
            });
        assert(allSame &&
               "All features must have the same referenceIDIndex to be coalesced into a "
               "GenomicFeatureGroup.");
        return refIDIndex;
    }
};

}  // namespace dataTypes
