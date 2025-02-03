// NOLINTBEGIN
#pragma once

#include <cassert>
#include <cstddef>
#include <unordered_set>
#include <vector>

#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "IITree.hpp"

namespace annotation {

using namespace dataTypes;

/**
 * @brief Parameters for feature merging operations.
 *
 * This structure holds configuration parameters for merging genomic features.
 * Merging operation of OPPOSITE is not allowed because it is a logical error.
 *
 * @param mergingOrientation The orientation to consider when merging features.
 * @param mergingTolerance The tolerance value defining proximity for merging.
 */
struct FeatureMergingParameters {
    FeatureMergingParameters(GenomicStrandSpecificity mergingSpecificity, int mergingTolerance)
        : mergingSpecificity(mergingSpecificity), mergingTolerance(mergingTolerance) {}

    GenomicStrandSpecificity mergingSpecificity;
    int mergingTolerance;
};

class GenomicFeatureTreeMerger {
   public:
    GenomicFeatureTreeMerger(const FeatureMergingParameters &parameters) : parameters(parameters) {}

    /**
     * @brief Merges overlapping genomic features within the provided feature tree.
     *
     * This method iterates through the genomic features in the feature tree,
     * identifies and merges features that are within the specified merging tolerance
     * and orientation. Invalidated feature indices are tracked and removed from the
     * feature tree after processing. Finally, the feature tree is re-indexed.
     *
     * @param featureTree The IITree containing GenomicFeature objects to be merged.
     */
    void merge(IITree<int, GenomicFeature> &featureTree);

   private:
    std::unordered_set<size_t> invalidatedIndices;
    FeatureMergingParameters parameters;

    static void updateIntervalsToMergedCoordinates(
        IITree<int, GenomicFeature> &featureTree) noexcept;

    [[nodiscard]] auto featureMeetsMergingConditions(
        const GenomicFeature &feature, const GenomicRegion &searchRegion) const noexcept -> bool;

    [[nodiscard]] auto updateSearchRegion(GenomicRegion &searchRegion,
                                          std::vector<size_t> &elementIndices,
                                          IITree<int, GenomicFeature> &featureTree) const noexcept
        -> bool;

    [[nodiscard]] auto mergeFeaturesWithinGenomicRegion(
        size_t baseFeatureIndex, IITree<int, GenomicFeature> &featureTree) const noexcept
        -> std::unordered_set<size_t>;
};

}  // namespace annotation

// NOLINTEND
