#include "GenomicFeatureTreeMerger.hpp"

// Standard
#include <algorithm>
#include <cstddef>
#include <unordered_set>
#include <utility>
#include <vector>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "GenomicStrandSpecificity.hpp"

namespace annotation {

void GenomicFeatureTreeMerger::merge(IITree<int, GenomicFeature> &featureTree) {
    if (featureTree.size() == 0) {
        return;
    }

    invalidatedIndices.clear();
    invalidatedIndices.reserve(featureTree.size());

    for (size_t index = 0; index < featureTree.size(); ++index) {
        if (invalidatedIndices.contains(index)) {
            continue;
        }

        auto mergingResults = mergeFeaturesWithinGenomicRegion(index, featureTree);
        invalidatedIndices.insert(mergingResults.begin(), mergingResults.end());
    }

    featureTree.remove(invalidatedIndices);

    updateIntervalsToMergedCoordinates(featureTree);

    featureTree.indexNoSort();
}

void GenomicFeatureTreeMerger::updateIntervalsToMergedCoordinates(
    IITree<int, GenomicFeature> &featureTree) noexcept {
    for (size_t index = 0; index < featureTree.size(); ++index) {
        const GenomicRegion &region = featureTree.getData(index).genomicRegion;

        featureTree.setIntervalStart(index, region.getStart());
        featureTree.setIntervalEnd(index, region.getEnd());
    }
}

auto GenomicFeatureTreeMerger::featureMeetsMergingConditions(
    const GenomicFeature &feature, const GenomicRegion &searchRegion) const noexcept -> bool {
    return searchRegion.overlapsWithTolerance(
        feature.genomicRegion,
        GenomicOrientation::fromStrandSpecificity(parameters.mergingSpecificity),
        parameters.mergingTolerance);
}

auto GenomicFeatureTreeMerger::updateSearchRegion(
    GenomicRegion &searchRegion, std::vector<size_t> &elementIndices,
    IITree<int, GenomicFeature> &featureTree) const noexcept -> bool {
    bool updatedRegion = false;
    std::vector<size_t> validatedFeatureIndices;
    validatedFeatureIndices.reserve(elementIndices.size());

    for (const size_t index : elementIndices) {
        const GenomicFeature &feature = featureTree.getData(index);

        if (featureMeetsMergingConditions(feature, searchRegion)) {
            validatedFeatureIndices.push_back(index);
            updatedRegion = searchRegion.combine(feature.genomicRegion) || updatedRegion;
        }
    }

    elementIndices = std::move(validatedFeatureIndices);
    return updatedRegion;
}

auto GenomicFeatureTreeMerger::mergeFeaturesWithinGenomicRegion(
    size_t baseFeatureIndex, IITree<int, GenomicFeature> &featureTree) const noexcept
    -> std::unordered_set<size_t> {
    GenomicFeature &baseFeature = featureTree.getData(baseFeatureIndex);
    GenomicRegion expandingSearchRegion = baseFeature.genomicRegion;

    if (parameters.mergingSpecificity == GenomicStrandSpecificity::UNSPECIFIC) {
        expandingSearchRegion.setStrand(GenomicStrand::NONE);
    }

    std::vector<size_t> overlappingFeatureIndices;
    bool updatedSearchRegion = false;

    do {
        overlappingFeatureIndices.clear();

        featureTree.overlap(expandingSearchRegion.getStart() - parameters.mergingTolerance,
                            expandingSearchRegion.getEnd() + parameters.mergingTolerance,
                            overlappingFeatureIndices);

        // Only contains original fragment
        if (overlappingFeatureIndices.size() <= 1) {
            return {};
        }

        // TODO: Check whether returned indices are always already sorted
        std::ranges::sort(overlappingFeatureIndices);

        // Remove the base element index since this is kept and updated to merge all
        // overlapping elements
        auto eraseIndices = std::ranges::remove(overlappingFeatureIndices, baseFeatureIndex);
        overlappingFeatureIndices.erase(eraseIndices.begin(), eraseIndices.end());

        updatedSearchRegion =
            updateSearchRegion(expandingSearchRegion, overlappingFeatureIndices, featureTree);
    } while (updatedSearchRegion);

    // Only set the internal start and end coordinate since setting the IITree
    // start end ends might cause issues with the IITree index
    baseFeature.genomicRegion.setStart(expandingSearchRegion.getStart());
    baseFeature.genomicRegion.setEnd(expandingSearchRegion.getEnd());
    baseFeature.genomicRegion.setStrand(expandingSearchRegion.getStrand());

    return {overlappingFeatureIndices.begin(), overlappingFeatureIndices.end()};
}

}  // namespace annotation
