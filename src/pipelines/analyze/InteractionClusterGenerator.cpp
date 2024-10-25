#include "InteractionClusterGenerator.hpp"

#include <algorithm>
#include <execution>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "GenomicRegion.hpp"
#include "InteractionSegment.hpp"
#include "Logger.hpp"

namespace pipelines::analyze {

/**
 * @brief Merges overlapping interaction clusters and returns these.
 *
 * This function takes a list of interaction clusters, sorts them from back to front,
 * and merges any overlapping clusters. This is done from back to front while closing clusters that
 * are further back than the current cluster.
 *
 * @param clusters A reference to the list of interaction clusters to be merged.
 * @return A list of finalized interaction clusters.
 */
auto InteractionClusterGenerator::mergeClusters(std::vector<InteractionCluster>&& clusters)
    -> InteractionClusterGenerator::Result {
    // Clusters should be sorted from back to front

    std::vector<InteractionCluster> localClusters = std::move(clusters);

    Logger::log(LogLevel::INFO, "(", sampleName, ") Sorting clusters");

    std::ranges::sort(localClusters, std::less<>{});

    Logger::log(LogLevel::INFO, "(", sampleName, ") Finished sorting clusters");

    for (auto& cluster : localClusters | std::views::reverse) {
        if (openClusterQueue.empty()) {
            openClusterQueue.emplace_front(std::move(cluster));
            continue;
        }

        bool clusterWasMerged = false;

        auto prevIter = openClusterQueue.before_begin();

        for (auto iter = openClusterQueue.begin(); iter != openClusterQueue.end();) {
            if (iter->overlaps(cluster, graceDistance)) {
                iter->merge(cluster);
                clusterWasMerged = true;
                break;
            }

            if (clusterIsBeforeOpenCluster(cluster, *iter)) {
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

    annotateSupplementaryFeatures();

    const size_t totalClusterCount = includedClusterCount + excludedClusterCount;
    Logger::log(LogLevel::INFO, "(", sampleName, ") Finished processing ", totalClusterCount,
                " clusters. Included ", includedClusterCount, " clusters, excluded ",
                excludedClusterCount, " clusters");

    return {.annotatedClusters = std::move(finishedClusters),
            .featureCounts = std::move(featureCounts),
            .supplementaryFeatureAnnotator = std::move(supplementaryFeatureAnnotator)};
}

void InteractionClusterGenerator::finalizeCluster(InteractionCluster&& cluster) noexcept {
    logClusteringStatus();

    if (cluster.fragmentCount() < minReadCount) {
        excludedClusterCount++;
        return;
    }

    const auto annotatedCluster = annotateCluster(std::move(cluster));

    // Already count the features for the annotated clusters
    // Partially annotated clusters will be counted when they get supplementary anntated
    if (const auto* annotated = std::get_if<AnnotatedInteractionCluster>(&annotatedCluster)) {
        finishedClusters.emplace_back(*annotated);
        featureCounts[annotated->getFirstFeatureID()]++;
        featureCounts[annotated->getSecondFeatureID()]++;
    } else {
        partiallyAnnotatedClusters.emplace_back(
            std::get<PartiallyAnnotatedInteractionCluster>(annotatedCluster));
    }

    includedClusterCount++;
}

auto InteractionClusterGenerator::clusterIsBeforeOpenCluster(
    const InteractionCluster& cluster, const InteractionCluster& openCluster) noexcept -> bool {
    return (cluster.getSecondSegment().getReferenceIDIndex() <
            openCluster.getSecondSegment().getReferenceIDIndex()) ||
           (cluster.getSecondSegment().getEnd() < openCluster.getSecondSegment().getStart());
}

auto InteractionClusterGenerator::annotateCluster(InteractionCluster&& cluster) noexcept
    -> std::variant<AnnotatedInteractionCluster, PartiallyAnnotatedInteractionCluster> {
    const auto& firstSegment = cluster.getFirstSegment();
    const std::string firstReferenceID = referenceIDs[firstSegment.getReferenceIDIndex()];

    const auto firstFeature = featureAnnotator->getBestOverlappingFeature(
        {firstReferenceID, firstSegment.getStart(), firstSegment.getEnd(),
         firstSegment.getStrand()},
        featureOrientation);

    const auto& secondSegment = cluster.getSecondSegment();
    const std::string secondReferenceID = referenceIDs[secondSegment.getReferenceIDIndex()];

    const auto secondFeature = featureAnnotator->getBestOverlappingFeature(
        {secondReferenceID, secondSegment.getStart(), secondSegment.getEnd(),
         secondSegment.getStrand()},
        featureOrientation);

    if (firstFeature && secondFeature) {
        return AnnotatedInteractionCluster{
            std::move(cluster), firstFeature.value().groupID.value_or(firstFeature.value().id),
            secondFeature.value().groupID.value_or(secondFeature.value().id)};
    }

    // Non annotated segments are inserted into the supplementary feature annotator
    // They are later merged and than non annotated segments are annotated by the supplementary
    // annotator
    if (!firstFeature.has_value()) {
        auto genomicRegion = GenomicRegion{firstReferenceID, firstSegment.getStart(),
                                           firstSegment.getEnd(), firstSegment.getStrand()};

        supplementaryFeatureAnnotator.insert(genomicRegion);
    }

    if (!secondFeature.has_value()) {
        auto genomicRegion = GenomicRegion{secondReferenceID, secondSegment.getStart(),
                                           secondSegment.getEnd(), secondSegment.getStrand()};

        supplementaryFeatureAnnotator.insert(genomicRegion);
    }

    return PartiallyAnnotatedInteractionCluster{
        std::move(cluster),
        firstFeature.has_value()
            ? std::make_optional(firstFeature.value().groupID.value_or(firstFeature.value().id))
            : std::nullopt,
        secondFeature.has_value()
            ? std::make_optional(secondFeature.value().groupID.value_or(secondFeature.value().id))
            : std::nullopt};
}

void InteractionClusterGenerator::annotateSupplementaryFeatures() noexcept {
    supplementaryFeatureAnnotator.mergeIndexAllOverlappingFeatures(-graceDistance);

    auto getFeatureIDForSegment = [&](const InteractionSegment& segment) -> std::string {
        const std::string referenceID = referenceIDs[segment.getReferenceIDIndex()];
        const auto feature = supplementaryFeatureAnnotator.getBestOverlappingFeature(
            {referenceID, segment.getStart(), segment.getEnd(), segment.getStrand()},
            featureOrientation);

        assert(feature.has_value() && "All segments should be annotated");

        return feature.value().groupID.value_or(feature.value().id);
    };

    for (auto& cluster : partiallyAnnotatedClusters) {
        const auto firstFeatureID = cluster.getFirstFeatureID().has_value()
                                        ? cluster.getFirstFeatureID().value()
                                        : getFeatureIDForSegment(cluster.getFirstSegment());

        const auto secondFeatureID = cluster.getSecondFeatureID().has_value()
                                         ? cluster.getSecondFeatureID().value()
                                         : getFeatureIDForSegment(cluster.getSecondSegment());

        featureCounts[firstFeatureID]++;
        featureCounts[secondFeatureID]++;

        finishedClusters.emplace_back(std::move(cluster), firstFeatureID, secondFeatureID);
    }
}

void InteractionClusterGenerator::logClusteringStatus() const noexcept {
    constexpr size_t LOGGING_INTERVAL = 100000;
    const size_t totalClusterCount = includedClusterCount + excludedClusterCount;

    if (totalClusterCount % LOGGING_INTERVAL != 0 || finishedClusters.empty()) [[likely]] {
        return;
    }

    Logger::log(LogLevel::INFO, "(", sampleName, ") Processed ", totalClusterCount,
                " clusters. Included ", includedClusterCount, " clusters, excluded ",
                excludedClusterCount, " clusters");
}

}  // namespace pipelines::analyze
