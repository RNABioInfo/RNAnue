#include "SplicingEvaluationStep.hpp"

// Standard
#include <cctype>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "SortedGenomicRegionPair.hpp"

namespace pipelines::detect {

auto SplicingEvaluationStep::getGroupedFeaturesAtBoundingRegions(
    const BoundingRegions &boundingRegions) const -> FeaturePairs {
    std::unordered_map<std::string, GenomicFeature> featuresFirstRegionByGroupID;

    // Iterate all features that overlap the upstream bounding region and check whether their end is
    // within upstream bounding region.
    for (const GenomicFeature &feature : config.featureAnnotator->overlappingFeatureIt(
             boundingRegions.first, config.featureOrientation)) {
        if (not feature.getParentID() ||
            !boundingRegions.first.contains(feature.getGenomicRegion().getEnd() - 1)) {
            continue;
        }

        featuresFirstRegionByGroupID.emplace(*feature.getParentID(), feature);
    }

    if (featuresFirstRegionByGroupID.empty()) {
        return {};
    }

    FeaturePairs groupedFeaturePairs;

    // Iterate all features that overlap the second bounding regions and check whether their start
    // is within the second bounding region.
    for (const GenomicFeature &feature : config.featureAnnotator->overlappingFeatureIt(
             boundingRegions.second, config.featureOrientation)) {
        if (!feature.getParentID() ||
            !featuresFirstRegionByGroupID.contains(*feature.getParentID()) ||
            !boundingRegions.second.contains(feature.getGenomicRegion().getStart())) {
            continue;
        }

        auto partnerFeature = featuresFirstRegionByGroupID.at(*feature.getParentID());

        if (feature.getGenomicRegion().getStrand() !=
            partnerFeature.getGenomicRegion().getStrand()) {
            continue;
        }

        groupedFeaturePairs.emplace_back(partnerFeature, feature);
    }

    return groupedFeaturePairs;
}

auto SplicingEvaluationStep::getBoundingRegionsForRecords(const SortedGenomicRegionPair &regionPair,
                                                          const size_t &tolerance)
    -> BoundingRegions {
    auto firstBoundingRegion = regionPair.firstRegion;
    firstBoundingRegion.setStart(firstBoundingRegion.getEnd() - 1);
    auto secondBoundingRegion = regionPair.secondRegion;
    secondBoundingRegion.setEnd(secondBoundingRegion.getStart() + 1);

    return std::make_pair(firstBoundingRegion.expanded(tolerance),
                          secondBoundingRegion.expanded(tolerance));
}

auto SplicingEvaluationStep::groupedFeaturePairsEncloseAdditonalFeatureFromGroup(
    const FeaturePairs &featurePairs) const -> bool {
    for (const FeaturePair &featurePair : featurePairs) {
        const GenomicRegion enclosingRegion{
            featurePair.first.getGenomicRegion().getReferenceIDIndex(),
            {.startPosition = featurePair.first.getGenomicRegion().getEnd(),
             .endPosition = featurePair.second.getGenomicRegion().getStart()},
            featurePair.first.getGenomicRegion().getStrand()};

        for (const GenomicFeature &enclosedFeature : config.featureAnnotator->overlappingFeatureIt(
                 enclosingRegion, config.featureOrientation)) {
            if (enclosedFeature.getParentID() == featurePair.first.getParentID()) {
                return true;
            }
        }
    }

    return false;
};

}  // namespace pipelines::detect
