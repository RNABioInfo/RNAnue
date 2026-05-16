#include "InteractionClusterGenerator.hpp"

// Standard
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "ClusteringParameters.hpp"
#include "GenomicFeature.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "InteractionCluster.hpp"
#include "InteractionClusterComponentBuilder.hpp"
#include "PartiallyAnnotatedInteractionCluster.hpp"
#include "VariantOverload.hpp"

namespace pipelines::analyze {

auto InteractionClusterGenerator::mergeClusters(std::vector<InteractionCluster> &&clusters)
    -> Result {
    std::vector<InteractionCluster> localClusters = std::move(clusters);

    auto clusterGroups = InteractionClusterComponentBuilder::groupClusters(
        std::move(localClusters), parameters.clusterMergingStrandSpecificity);

    for (auto &clusterGroup : clusterGroups) {
        finalizeMergedClusters(
            InteractionClusterComponentBuilder::mergeGroup(std::move(clusterGroup), parameters));
    }

    return releaseResult();
};

auto InteractionClusterGenerator::mergeClusterGroup(std::vector<InteractionCluster> &&clusters)
    -> Result {
    finalizeMergedClusters(
        InteractionClusterComponentBuilder::mergeGroup(std::move(clusters), parameters));
    return releaseResult();
};

auto InteractionClusterGenerator::clustersOverlap(const InteractionCluster &cluster1,
                                                  const InteractionCluster &cluster2,
                                                  const ClusteringParameters &parameters) noexcept
    -> bool {
    return InteractionClusterComponentBuilder::clustersOverlap(cluster1, cluster2, parameters);
};

auto InteractionClusterGenerator::annotateCluster(InteractionCluster &&cluster) noexcept
    -> std::variant<AnnotatedInteractionCluster, PartiallyAnnotatedInteractionCluster> {
    const auto &firstSegment = cluster.getFirstSegment();
    const auto &secondSegment = cluster.getSecondSegment();

    const auto firstFeature =
        featureAnnotator->getBestOverlappingFeature(firstSegment, parameters.featureOrientation);

    const auto secondFeature =
        featureAnnotator->getBestOverlappingFeature(secondSegment, parameters.featureOrientation);

    auto extractFeatureID =
        [](const std::optional<GenomicFeature> &feature) -> std::optional<std::string> {
        if (!feature) {
            return std::nullopt;
        }
        // Prefer groupID if it exists; otherwise fallback to featureID
        return feature->getAnnotationID();
    };

    auto firstFeatureID = extractFeatureID(firstFeature);
    auto secondFeatureID = extractFeatureID(secondFeature);

    // If fully annotated, return the annotated cluster
    if (firstFeatureID && secondFeatureID) {
        return AnnotatedInteractionCluster{std::move(cluster), *firstFeatureID, *secondFeatureID};
    }

    // Otherwise, populate supplementary features if the cluster passes filters
    if (clusterPassesFilters(cluster)) {
        if (!firstFeatureID) {
            // Create a new supplementary feature from the first segment

            // If the feature orientation is opposite switch strand for segment
            GenomicRegion firstRegionCopy = cluster.getFirstSegment();
            if (parameters.featureOrientation == GenomicOrientation::OPPOSITE) {
                firstRegionCopy.setStrand(!firstRegionCopy.getStrand());
            }

            supplementaryFeatureRegions[firstSegment.getReferenceIDIndex()].emplace_back(
                asGenomicFeature(firstRegionCopy));
        }

        if (!secondFeatureID) {
            // Create a new supplementary feature from the second segment

            // If the feature orientation is opposite switch strand for segment
            GenomicRegion secondRegionCopy = cluster.getSecondSegment();
            if (parameters.featureOrientation == GenomicOrientation::OPPOSITE) {
                secondRegionCopy.setStrand(!secondRegionCopy.getStrand());
            }
            supplementaryFeatureRegions[secondSegment.getReferenceIDIndex()].emplace_back(
                asGenomicFeature(secondRegionCopy));
        }
    }

    // Return partially annotated cluster
    return PartiallyAnnotatedInteractionCluster{std::move(cluster), std::move(firstFeatureID),
                                                std::move(secondFeatureID)};
};

auto InteractionClusterGenerator::clusterPassesFilters(
    const InteractionCluster &cluster) const noexcept -> bool {
    return cluster.getTranscriptContribution() >= parameters.minimumClusterTrascriptContribution &&
           cluster.segmentsMaxSelfOverlapFraction() <= parameters.maxClusterSelfOverlapFraction;
};

void InteractionClusterGenerator::attributeCluster(AnnotatedInteractionCluster &&cluster) noexcept {
    // Update counts

    featureCountsByFeatureID[cluster.getFirstFeatureID()] += cluster.getTranscriptContribution();
    featureCountsByFeatureID[cluster.getSecondFeatureID()] += cluster.getTranscriptContribution();

    if (clusterPassesFilters(cluster)) {
        finishedClusters.emplace_back(std::move(cluster));
        ++includedClusterCount;
    } else {
        ++excludedClusterCount;
    }
}
void InteractionClusterGenerator::attributeCluster(
    PartiallyAnnotatedInteractionCluster &&cluster) noexcept {
    if (clusterPassesFilters(cluster)) {
        partiallyAnnotatedClusters.emplace_back(std::move(cluster));
        ++includedClusterCount;
    } else {
        ++excludedClusterCount;
    }
}
void InteractionClusterGenerator::finalizeCluster(InteractionCluster &&cluster) noexcept {
    // Perform annotation and attribution in one place
    std::visit(
        [&](auto &&resultingCluster) {
            attributeCluster(std::forward<decltype(resultingCluster)>(resultingCluster));
        },
        annotateCluster(std::forward<InteractionCluster>(cluster)));
};

void InteractionClusterGenerator::finalizeMergedClusters(
    std::vector<InteractionCluster> &&clusters) noexcept {
    for (auto &cluster : clusters) {
        finalizeCluster(std::move(cluster));
    }
}

auto InteractionClusterGenerator::releaseResult() noexcept -> Result {
    return {.finishedClusters = std::move(finishedClusters),
            .partiallyAnnotatedClusters = std::move(partiallyAnnotatedClusters),
            .supplementaryFeatureMap = std::move(supplementaryFeatureRegions),
            .featureCounts = std::move(featureCountsByFeatureID),
            .includedClusterCount = includedClusterCount,
            .excludedClusterCount = excludedClusterCount};
}

void InteractionClusterGenerator::Result::merge(Result &&other) noexcept {
    // Merge finished clusters
    finishedClusters.reserve(finishedClusters.size() + other.finishedClusters.size());
    finishedClusters.insert(finishedClusters.end(),
                            std::make_move_iterator(other.finishedClusters.begin()),
                            std::make_move_iterator(other.finishedClusters.end()));

    // Merge partially annotated clusters
    partiallyAnnotatedClusters.reserve(partiallyAnnotatedClusters.size() +
                                       other.partiallyAnnotatedClusters.size());
    partiallyAnnotatedClusters.insert(
        partiallyAnnotatedClusters.end(),
        std::make_move_iterator(other.partiallyAnnotatedClusters.begin()),
        std::make_move_iterator(other.partiallyAnnotatedClusters.end()));

    // Merge supplementary feature maps
    for (auto &[key, features] : other.supplementaryFeatureMap) {
        auto &thisFeatures = supplementaryFeatureMap[key];
        thisFeatures.reserve(thisFeatures.size() + features.size());
        thisFeatures.insert(thisFeatures.end(), std::make_move_iterator(features.begin()),
                            std::make_move_iterator(features.end()));
    }

    // Merge feature counts
    for (const auto &[featureID, count] : other.featureCounts) {
        featureCounts[featureID] += count;
    }

    // Finally update counters
    includedClusterCount += other.includedClusterCount;
    excludedClusterCount += other.excludedClusterCount;
}

}  // namespace pipelines::analyze
