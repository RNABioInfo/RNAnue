#include "ParallelInteractionClusterGenerator.hpp"

// Standard
#include <utils/strings.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iterator>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

// Internal
#include "FeatureAnnotator.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "InteractionCluster.hpp"
#include "InteractionClusterComponentBuilder.hpp"
#include "InteractionClusterGenerator.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"

namespace pipelines::analyze {

auto ParallelInteractionClusterGenerator::mergeClusters(std::vector<InteractionCluster>&& clusters,
                                                        size_t threadCount, size_t batchSize)
    -> Result {
    std::mutex mergeMutex;
    std::mutex groupMutex;
    (void)batchSize;

    std::vector<InteractionCluster> localClusters = std::move(clusters);

    Logger::log("Grouping clusters");

    auto clusterGroups = InteractionClusterComponentBuilder::groupClusters(
        std::move(localClusters), parameters.clusterMergingStrandSpecificity);
    size_t nextGroupIndex = 0;

    Logger::log("Finished grouping clusters into ", clusterGroups.size(), " independent groups");

    auto clusterGroupConsumer = [&]() {
        while (true) {
            std::vector<InteractionCluster> clusterGroup;

            {
                std::lock_guard<std::mutex> lock(groupMutex);
                if (nextGroupIndex >= clusterGroups.size()) {
                    break;
                }

                clusterGroup = std::move(clusterGroups[nextGroupIndex]);
                ++nextGroupIndex;
            }

            InteractionClusterGenerator clusterGenerator{featureAnnotator, parameters};

            auto result = clusterGenerator.mergeClusterGroup(std::move(clusterGroup));

            {
                std::lock_guard<std::mutex> lock(mergeMutex);
                clusteringResults.merge(std::move(result));
                logClusteringStatus();
            }
        }
    };

    std::vector<std::thread> consumerThreads;
    const size_t workerCount = (std::max)(size_t{1}, threadCount);
    consumerThreads.reserve(workerCount);

    for (size_t i = 0; i < workerCount; ++i) {
        consumerThreads.emplace_back(clusterGroupConsumer);
    }

    for (std::thread& consumerThread : consumerThreads) {
        consumerThread.join();
    }

    FeatureAnnotator supplementaryFeatureAnnotator{clusteringResults.supplementaryFeatureMap};

    supplementaryFeatureAnnotator.mergeAllOverlappingFeatures(
        {parameters.featureOrientation.strandSpecificity(), 1});

    annotatePartiallyAnnotatedClusters(supplementaryFeatureAnnotator);

    Logger::log("Finished clustering. Total clusters after cluster contributions score filtering: ",
                clusteringResults.finishedClusters.size());

    return {.annotatedClusters = std::move(clusteringResults.finishedClusters),
            .featureCounts = std::move(clusteringResults.featureCounts),
            .supplementaryFeatureAnnotator = std::move(supplementaryFeatureAnnotator)};
}

void ParallelInteractionClusterGenerator::annotatePartiallyAnnotatedClusters(
    const FeatureAnnotator& supplementaryFeatureAnnotator) noexcept {
    auto getSupplementaryFeatureIDForSegment = [&](const GenomicRegion& segment) -> std::string {
        auto features = supplementaryFeatureAnnotator.getOverlappingFeatures(
            segment, parameters.featureOrientation);

        // If the number of features is not exactly one and merging is unspecific,
        // attempt to narrow down the features by selecting the best overlapping one.
        if (features.size() != 1 &&
            parameters.clusterMergingStrandSpecificity == GenomicStrandSpecificity::UNSPECIFIC) {
            if (auto bestFeature =
                    supplementaryFeatureAnnotator.getBestOverlappingFeatureWithPreferredOrientation(
                        segment, parameters.featureOrientation)) {
                features = {*bestFeature};  // Use the best match as the unique feature.
            }
        }

        // Warn if we still don't have exactly one feature.
        if (features.size() != 1) {
            Logger::log<LogLevel::WARNING>(
                "Expected exactly one supplementary feature hit but got: ");
            for (const auto& feature : features) {
                Logger::log<LogLevel::WARNING>(feature);
            }

            auto otherFeatures = supplementaryFeatureAnnotator.getOverlappingFeatures(
                segment, GenomicOrientation::BOTH);

            Logger::log<LogLevel::WARNING>("Region is: ", segment,
                                           ", Other features in region are: ");
            for (const auto& feature : otherFeatures) {
                Logger::log<LogLevel::WARNING>(feature);
            }
        }

        // Ensure the invariant with an assert.
        assert(features.size() == 1 &&
               "All segments should be annotated and have a unique feature associated");

        return features.front().getAnnotationID();
    };

    for (auto& partiallyAnnotatedCluster : clusteringResults.partiallyAnnotatedClusters) {
        std::string firstFeatureID;
        if (partiallyAnnotatedCluster.getFirstFeatureID().has_value()) {
            firstFeatureID = partiallyAnnotatedCluster.getFirstFeatureID().value();
        } else {
            firstFeatureID =
                getSupplementaryFeatureIDForSegment(partiallyAnnotatedCluster.getFirstSegment());
        }

        std::string secondFeatureID;
        if (partiallyAnnotatedCluster.getSecondFeatureID().has_value()) {
            secondFeatureID = partiallyAnnotatedCluster.getSecondFeatureID().value();
        } else {
            secondFeatureID =
                getSupplementaryFeatureIDForSegment(partiallyAnnotatedCluster.getSecondSegment());
        }

        clusteringResults.featureCounts[firstFeatureID] +=
            partiallyAnnotatedCluster.getTranscriptContribution();
        clusteringResults.featureCounts[secondFeatureID] +=
            partiallyAnnotatedCluster.getTranscriptContribution();

        if (partiallyAnnotatedCluster.getTranscriptContribution() >=
                parameters.minimumClusterTrascriptContribution &&
            partiallyAnnotatedCluster.segmentsMaxSelfOverlapFraction() <=
                parameters.maxClusterSelfOverlapFraction) {
            clusteringResults.finishedClusters.emplace_back(partiallyAnnotatedCluster,
                                                            firstFeatureID, secondFeatureID);
        }
    }
}

void ParallelInteractionClusterGenerator::logClusteringStatus() const noexcept {
    Logger::log("Processed ", clusteringResults.totalClusterCount(), " clusters. Included ",
                clusteringResults.includedClusterCount, " clusters, excluded ",
                clusteringResults.excludedClusterCount, " clusters");
}

}  // namespace pipelines::analyze
