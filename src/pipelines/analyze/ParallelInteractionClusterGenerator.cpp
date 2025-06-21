#include "ParallelInteractionClusterGenerator.hpp"

// Standard
#include <utils/strings.h>

#include <algorithm>
#include <cassert>
#include <condition_variable>
#include <cstddef>
#include <functional>
#include <iterator>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <utility>
#include <vector>

// Internal
#include "FeatureAnnotator.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "InteractionCluster.hpp"
#include "InteractionClusterGenerator.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"

namespace pipelines::analyze {

auto ParallelInteractionClusterGenerator::mergeClusters(std::vector<InteractionCluster>&& clusters,
                                                        size_t threadCount, size_t batchSize)
    -> Result {
    std::mutex produceMutex;
    std::mutex mergeMutex;
    std::condition_variable conditionVariable;
    bool doneProducing = false;

    std::vector<InteractionCluster> localClusters = std::move(clusters);

    Logger::log("Sorting clusters");

    // Clusters should be sorted from back to front
    std::ranges::sort(localClusters, std::less<>{});

    Logger::log("Finished sorting clusters");

    std::queue<std::vector<InteractionCluster>> clusterBatches;

    auto clusterBatchProducer = [&]() {
        while (!localClusters.empty()) {
            auto startIterator = localClusters.size() > batchSize
                                     ? localClusters.end() - long(batchSize)
                                     : localClusters.begin();

            // Check that the first cluster in the batch does not overlap with the last cluster
            // in the remaining clusters
            while (startIterator != localClusters.begin() &&
                   !(startIterator - 1)->isBefore(*startIterator)) {
                startIterator--;
            }

            const size_t batchElementCount = std::distance(startIterator, localClusters.end());
            auto batch = std::vector<InteractionCluster>();
            batch.reserve(batchElementCount);

            auto endIterator = localClusters.end();

            // Move elements from the localClusters vector to the batch vector and erase them
            // from the localClusters vector
            std::move(startIterator, endIterator, std::back_inserter(batch));
            localClusters.erase(startIterator, endIterator);

            {
                std::unique_lock<std::mutex> lock(produceMutex);
                clusterBatches.push(std::move(batch));
                conditionVariable.notify_all();
            }
        }

        {
            std::unique_lock<std::mutex> lock(produceMutex);
            doneProducing = true;
            conditionVariable.notify_all();
        }
    };

    auto clusterBatchConsumer = [&]() {
        while (true) {
            std::vector<InteractionCluster> batch;

            {
                std::unique_lock<std::mutex> lock(produceMutex);
                conditionVariable.wait(lock,
                                       [&] { return !clusterBatches.empty() || doneProducing; });

                if (clusterBatches.empty() && doneProducing) {
                    break;
                }

                batch = std::move(clusterBatches.front());
                clusterBatches.pop();
            }

            InteractionClusterGenerator clusterGenerator{featureAnnotator, parameters};

            auto result = clusterGenerator.mergeClusters(std::move(batch));

            {
                std::lock_guard<std::mutex> lock(mergeMutex);
                clusteringResults.merge(std::move(result));
                logClusteringStatus();
            }
        }
    };

    std::thread producerThread(clusterBatchProducer);

    std::vector<std::thread> consumerThreads;
    consumerThreads.reserve(threadCount);

    for (size_t i = 0; i < threadCount; ++i) {
        consumerThreads.emplace_back(clusterBatchConsumer);
    }

    producerThread.join();

    for (std::thread& consumerThread : consumerThreads) {
        consumerThread.join();
    }

    FeatureAnnotator supplementaryFeatureAnnotator{clusteringResults.supplementaryFeatureMap};

    supplementaryFeatureAnnotator.mergeAllOverlappingFeatures(
        {parameters.featureOrientation.strandSpecificity(), 1});

    annotatePartiallyAnnotatedClusters(supplementaryFeatureAnnotator);

    Logger::log("Finished clustering");

    return {.annotatedClusters = std::move(clusteringResults.finishedClusters),
            .featureCounts = std::move(clusteringResults.featureCounts),
            .supplementaryFeatureAnnotator = std::move(supplementaryFeatureAnnotator)};
}

void ParallelInteractionClusterGenerator::annotatePartiallyAnnotatedClusters(
    const FeatureAnnotator& supplementaryFeatureAnnotator) noexcept {
    auto getSupplementaryFeatureIDForSegment = [&](const GenomicRegion& segment) -> std::string {
        const auto features = supplementaryFeatureAnnotator.getOverlappingFeatures(
            segment,
            GenomicOrientation::fromStrandSpecificity(parameters.clusterMergingStrandSpecificity));

        // TODO: Fix this sometimes not working with strand specific
        if (features.empty()) {
            Logger::log<LogLevel::WARNING>("Could not find supplementary annotation for region: ",
                                           segment, "Found following features: ");

            const auto features = supplementaryFeatureAnnotator.getOverlappingFeatures(
                segment, GenomicOrientation::BOTH);

            for (const auto& feature : features) {
                Logger::log<LogLevel::WARNING>(feature);
            }
        }
        assert((features.size() == 1) &&
               "All segments should be annotated and should have a unique feature associated");

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
