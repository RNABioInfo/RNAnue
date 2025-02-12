#include "SplitRecordsSplicingEvaluator.hpp"

// Standard
#include <cctype>
#include <cstddef>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "SortedGenomicRegionPair.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"

auto SplitRecordsSplicingEvaluator::isSplicedSplitRecord(
    const SplitRecords &splitRecords,
    const SplitRecordsEvaluationParameters::SplicingParameters &parameters) -> bool {
    const std::optional<GenomicRegion> regionOne = GenomicRegion::fromSamRecord(splitRecords[0]);
    const std::optional<GenomicRegion> regionTwo = GenomicRegion::fromSamRecord(splitRecords[1]);

    if (!regionOne.has_value() || !regionTwo.has_value()) {
        return false;
    }

    const auto boundingRegions =
        getBoundingRegionsForRecords({*regionOne, *regionTwo}, parameters.splicingTolerance);

    const auto featurePairs = getGroupedFeaturesAtBoundingRegions(boundingRegions, parameters);

    if (parameters.allowAlternativeSplicing || featurePairs.empty()) {
        return !featurePairs.empty();
    }

    return !groupedFeaturePairsEncloseAdditonalFeatureFromGroup(featurePairs, parameters);
}

auto SplitRecordsSplicingEvaluator::getGroupedFeaturesAtBoundingRegions(
    const BoundingRegions &boundingRegions,
    const SplitRecordsEvaluationParameters::SplicingParameters &parameters) -> FeaturePairs {
    std::unordered_map<std::string, GenomicFeature> featuresFirstRegionByGroupID;

    // Iterate all features that overlap the upstream bounding region and check whether their end is
    // within upstream bounding region.
    for (const GenomicFeature &feature : parameters.featureAnnotator->overlappingFeatureIt(
             boundingRegions.first, parameters.orientation)) {
        if (not feature.getGroupID() ||
            !boundingRegions.first.contains(feature.getGenomicRegion().getEnd() - 1)) {
            continue;
        }

        featuresFirstRegionByGroupID.emplace(*feature.getGroupID(), feature);
    }

    if (featuresFirstRegionByGroupID.empty()) {
        return {};
    }

    FeaturePairs groupedFeaturePairs;

    // Iterate all features that overlap the second bounding regions and check whether their start
    // is within the second bounding region.
    for (const GenomicFeature &feature : parameters.featureAnnotator->overlappingFeatureIt(
             boundingRegions.second, parameters.orientation)) {
        if (!feature.getGroupID() ||
            !featuresFirstRegionByGroupID.contains(*feature.getGroupID()) ||
            !boundingRegions.second.contains(feature.getGenomicRegion().getStart())) {
            continue;
        }

        auto partnerFeature = featuresFirstRegionByGroupID.at(*feature.getGroupID());

        if (feature.getGenomicRegion().getStrand() !=
            partnerFeature.getGenomicRegion().getStrand()) {
            continue;
        }

        groupedFeaturePairs.emplace_back(partnerFeature, feature);
    }

    return groupedFeaturePairs;
}

auto SplitRecordsSplicingEvaluator::getBoundingRegionsForRecords(
    const SortedGenomicRegionPair &regionPair, const size_t &tolerance) -> BoundingRegions {
    auto firstBoundingRegion = regionPair.firstRegion;
    firstBoundingRegion.setStart(firstBoundingRegion.getEnd() - 1);
    auto secondBoundingRegion = regionPair.secondRegion;
    secondBoundingRegion.setEnd(secondBoundingRegion.getStart() + 1);

    return std::make_pair(firstBoundingRegion.expanded(tolerance),
                          secondBoundingRegion.expanded(tolerance));
}

auto SplitRecordsSplicingEvaluator::groupedFeaturePairsEncloseAdditonalFeatureFromGroup(
    const FeaturePairs &featurePairs,
    const SplitRecordsEvaluationParameters::SplicingParameters &parameters) -> bool {
    for (const FeaturePair &featurePair : featurePairs) {
        const GenomicRegion enclosingRegion{
            featurePair.first.getGenomicRegion().getReferenceIDIndex(),
            {.startPosition = featurePair.first.getGenomicRegion().getEnd(),
             .endPosition = featurePair.second.getGenomicRegion().getStart()},
            featurePair.first.getGenomicRegion().getStrand()};

        for (const GenomicFeature &enclosedFeature :
             parameters.featureAnnotator->overlappingFeatureIt(enclosingRegion,
                                                               parameters.orientation)) {
            if (enclosedFeature.getGroupID() == featurePair.first.getGroupID()) {
                return true;
            }
        }
    }

    return false;
};
