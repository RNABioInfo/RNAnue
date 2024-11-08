#include "InteractionClusterGenerator.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "GenomicRegion.hpp"
#include "InteractionSegment.hpp"

namespace pipelines::analyze {

/**
 * @brief Merges overlapping interaction clusters and returns these.
 *
 * This function takes a sorted list of interaction clusters,
 * and merges any overlapping clusters. This is done from back to front while closing clusters that
 * are further back than the current cluster.
 *
 * @param clusters A forwaring reference to the sorted list of interaction clusters to be merged.
 * @return A list of finalized interaction clusters.
 */
auto InteractionClusterGenerator::mergeClusters(std::vector<InteractionCluster>&& clusters)
    -> InteractionClusterGenerator::Result {
    std::vector<InteractionCluster> localClusters = std::move(clusters);

    for (auto& cluster : localClusters | std::views::reverse) {
        if (openClusterQueue.empty()) {
            openClusterQueue.emplace_front(std::move(cluster));
            continue;
        }

        bool clusterWasMerged = false;

        auto prevIter = openClusterQueue.before_begin();

        for (auto iter = openClusterQueue.begin(); iter != openClusterQueue.end();) {
            if (iter->overlaps(cluster, parameters.graceDistance)) {
                iter->merge(cluster);
                clusterWasMerged = true;
                break;
            }

            if (cluster.isBefore(*iter)) {
                finalizeCluster(std::move(*iter));
                iter = openClusterQueue.erase_after(prevIter);
                continue;
            }

            prevIter = iter;
            ++iter;
        }

        if (!clusterWasMerged) {
            openClusterQueue.emplace_after(prevIter, std::move(cluster));
        }
    }

    for (auto& cluster : openClusterQueue) {
        finalizeCluster(std::move(cluster));
    }

    return {.finishedClusters = std::move(finishedClusters),
            .partiallyAnnotatedClusters = std::move(partiallyAnnotatedClusters),
            .supplementaryFeatureRegions = std::move(supplementaryFeatureRegions),
            .featureCounts = std::move(featureCounts),
            .includedClusterCount = includedClusterCount,
            .excludedClusterCount = excludedClusterCount};
}

void InteractionClusterGenerator::finalizeCluster(InteractionCluster&& cluster) noexcept {
    auto annotatedCluster = annotateCluster(std::move(cluster));

    std::visit(
        [this](auto&& cluster) {
            using T = std::decay_t<decltype(cluster)>;
            if constexpr (std::is_same_v<T, AnnotatedInteractionCluster>) {
                featureCounts[cluster.getFirstFeatureID()] += cluster.fragmentCount();
                featureCounts[cluster.getSecondFeatureID()] += cluster.fragmentCount();

                if (cluster.fragmentCount() >= parameters.minReadCount &&
                    cluster.segmentsMaxOverlapFraction() <= parameters.maxOverlapFraction) {
                    finishedClusters.emplace_back(std::forward<decltype(cluster)>(cluster));
                    ++includedClusterCount;
                } else {
                    ++excludedClusterCount;
                }
            } else if constexpr (std::is_same_v<T, PartiallyAnnotatedInteractionCluster>) {
                if (cluster.fragmentCount() >= parameters.minReadCount &&
                    cluster.segmentsMaxOverlapFraction() <= parameters.maxOverlapFraction) {
                    ++includedClusterCount;
                } else {
                    ++excludedClusterCount;
                }
                partiallyAnnotatedClusters.emplace_back(std::forward<decltype(cluster)>(cluster));
            }
        },
        std::move(annotatedCluster));
}

auto InteractionClusterGenerator::annotateCluster(InteractionCluster&& cluster) noexcept
    -> std::variant<AnnotatedInteractionCluster, PartiallyAnnotatedInteractionCluster> {
    const auto& firstSegment = cluster.getFirstSegment();
    const auto& secondSegment = cluster.getSecondSegment();

    const auto firstReferenceID = referenceIDs[firstSegment.getReferenceIDIndex()];
    const auto secondReferenceID = referenceIDs[secondSegment.getReferenceIDIndex()];

    const auto firstFeature = featureAnnotator->getBestOverlappingFeature(
        {firstReferenceID, firstSegment.getStart(), firstSegment.getEnd(),
         firstSegment.getStrand()},
        parameters.featureOrientation);

    const auto secondFeature = featureAnnotator->getBestOverlappingFeature(
        {secondReferenceID, secondSegment.getStart(), secondSegment.getEnd(),
         secondSegment.getStrand()},
        parameters.featureOrientation);

    auto getFeatureID = [](const auto& feature) -> std::optional<std::string> {
        return feature ? std::make_optional(feature->groupID.value_or(feature->id)) : std::nullopt;
    };

    auto firstFeatureID = getFeatureID(firstFeature);
    auto secondFeatureID = getFeatureID(secondFeature);

    if (firstFeatureID && secondFeatureID) {
        return AnnotatedInteractionCluster{std::move(cluster), *firstFeatureID, *secondFeatureID};
    }

    if (!firstFeatureID) {
        supplementaryFeatureRegions.emplace_back(firstReferenceID, firstSegment.getStart(),
                                                 firstSegment.getEnd(), firstSegment.getStrand());
    }

    if (!secondFeatureID) {
        supplementaryFeatureRegions.emplace_back(secondReferenceID, secondSegment.getStart(),
                                                 secondSegment.getEnd(), secondSegment.getStrand());
    }

    return PartiallyAnnotatedInteractionCluster{std::move(cluster), std::move(firstFeatureID),
                                                std::move(secondFeatureID)};
}

}  // namespace pipelines::analyze
