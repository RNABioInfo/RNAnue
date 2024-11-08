#include "ParallelInteractionClusterGenerator.hpp"

// Standard
#include <algorithm>
#include <condition_variable>
#include <cstddef>
#include <iterator>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

// Internal
#include "GenomicRegion.hpp"
#include "InteractionCluster.hpp"
#include "InteractionClusterGenerator.hpp"
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
auto ParallelInteractionClusterGenerator::mergeClusters(std::vector<InteractionCluster>&& clusters,
                                                        const size_t threadCount,
                                                        const size_t batchSize)
    -> ParallelInteractionClusterGenerator::Result {
    std::mutex produceMutex;
    std::mutex mergeMutex;
    std::condition_variable conditionVariable;
    bool doneProducing = false;

    std::vector<InteractionCluster> localClusters = std::move(clusters);

    Logger::log(LogLevel::INFO, "Sorting clusters");

    // Clusters should be sorted from back to front
    std::ranges::sort(localClusters, std::less<>{});

    Logger::log(LogLevel::INFO, "Finished sorting clusters");

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

            InteractionClusterGenerator clusterGenerator{featureAnnotator, referenceIDs,
                                                         parameters};

            auto result = clusterGenerator.mergeClusters(std::move(batch));

            {
                std::lock_guard<std::mutex> lock(mergeMutex);
                mergeClusteringResults(std::move(result));
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

    annotateSupplementaryFeatures();

    const size_t totalClusterCount = includedClusterCount + excludedClusterCount;
    Logger::log(LogLevel::INFO, "Finished processing ", totalClusterCount, " clusters. Included ",
                includedClusterCount, " clusters, excluded ", excludedClusterCount, " clusters");

    return {.annotatedClusters = std::move(finishedClusters),
            .featureCounts = std::move(featureCounts),
            .supplementaryFeatureAnnotator = std::move(supplementaryFeatureAnnotator)};
}

void ParallelInteractionClusterGenerator::mergeClusteringResults(
    InteractionClusterGenerator::Result&& result) noexcept {
    InteractionClusterGenerator::Result localResult = std::move(result);

    finishedClusters.reserve(finishedClusters.size() + localResult.finishedClusters.size());
    finishedClusters.insert(finishedClusters.end(),
                            std::make_move_iterator(localResult.finishedClusters.begin()),
                            std::make_move_iterator(localResult.finishedClusters.end()));

    partiallyAnnotatedClusters.reserve(partiallyAnnotatedClusters.size() +
                                       localResult.partiallyAnnotatedClusters.size());
    partiallyAnnotatedClusters.insert(
        partiallyAnnotatedClusters.end(),
        std::make_move_iterator(localResult.partiallyAnnotatedClusters.begin()),
        std::make_move_iterator(localResult.partiallyAnnotatedClusters.end()));

    for (const GenomicRegion& region : localResult.supplementaryFeatureRegions) {
        supplementaryFeatureAnnotator.insert(region);
    }

    for (const auto& [featureID, count] : localResult.featureCounts) {
        featureCounts[featureID] += count;
    }

    includedClusterCount += localResult.includedClusterCount;
    excludedClusterCount += localResult.excludedClusterCount;
}

void ParallelInteractionClusterGenerator::annotateSupplementaryFeatures() noexcept {
    supplementaryFeatureAnnotator.mergeIndexAllOverlappingFeatures(-parameters.graceDistance);

    auto getFeatureIDForSegment = [&](const InteractionSegment& segment) -> std::string {
        const std::string referenceID = referenceIDs[segment.getReferenceIDIndex()];
        const auto feature = supplementaryFeatureAnnotator.getBestOverlappingFeature(
            {referenceID, segment.getStart(), segment.getEnd(), segment.getStrand()},
            parameters.featureOrientation);

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

        if (cluster.fragmentCount() >= parameters.minReadCount &&
            cluster.segmentsMaxOverlapFraction() <= parameters.maxOverlapFraction) {
            finishedClusters.emplace_back(std::move(cluster), firstFeatureID, secondFeatureID);
        }
    }
}

void ParallelInteractionClusterGenerator::logClusteringStatus() const noexcept {
    const size_t totalClusterCount = includedClusterCount + excludedClusterCount;
    Logger::log(LogLevel::INFO, "Processed ", totalClusterCount, " clusters. Included ",
                includedClusterCount, " clusters, excluded ", excludedClusterCount, " clusters");
}

}  // namespace pipelines::analyze
