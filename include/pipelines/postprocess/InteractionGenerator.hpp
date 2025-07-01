#pragma once

// Standard

// Internal
#include <algorithm>
#include <functional>
#include <list>
#include <ranges>
#include <utility>
#include <vector>

#include "GenomicStrandSpecificity.hpp"
#include "Interaction.hpp"
#include "Logger.hpp"

namespace pipelines::postprocess {

class InteractionGenerator {
   public:
    struct ClusteringParameters {
        float minSegmentOverlapFraction;
        GenomicStrandSpecificity mergingStrandSpecificity;
    };

    struct Result {
        std::vector<Interaction> superInteractions;
    };

    explicit InteractionGenerator(ClusteringParameters parameters) : parameters(parameters) {}

    [[nodiscard]] auto merge(std::vector<Interaction> &&interactions) noexcept -> Result {
        std::vector<Interaction> localClusters{std::move(interactions)};

        Logger::log("Sorting Interactions");
        std::ranges::sort(localClusters, std::less<>{});
        Logger::log("Finished sorting interactions");

        // Process from the end to the beginning
        for (auto &cluster : localClusters | std::views::reverse) {
            if (openClusterQueue.empty()) {
                openClusterQueue.emplace_front(std::move(cluster));
                continue;
            }

            // Newest leftmost clusters shall always be at the front
            if (cluster.isBefore(openClusterQueue.front())) {
                for (auto &cluster : openClusterQueue) {
                    finishedClusters.emplace_back(std::move(cluster));
                }

                openClusterQueue.clear();
            }

            bool clusterMerged = false;

            // Try merging with clusters in the open queue
            for (auto iter = openClusterQueue.begin(); iter != openClusterQueue.end();) {
                // TODO: Fix bug where multiple overlaps result not in merging into one cluster

                if (clustersOverlap(*iter, cluster, parameters) &&
                    iter->merge(cluster, parameters.mergingStrandSpecificity)) {
                    clusterMerged = true;
                    greedyMerge(iter);
                    break;
                }

                ++iter;
            }

            // If not merged, insert new cluster into the queue
            if (!clusterMerged) {
                openClusterQueue.emplace_front(std::move(cluster));
            }
        }

        // Finalize remaining clusters
        for (auto &cluster : openClusterQueue) {
            finishedClusters.emplace_back(std::move(cluster));
        }

        openClusterQueue.clear();

        return {.superInteractions = std::move(finishedClusters)};
    }

   private:
    ClusteringParameters parameters;
    std::list<Interaction> openClusterQueue;
    std::vector<Interaction> finishedClusters;

    static auto clustersOverlap(const Interaction &cluster1, const Interaction &cluster2,
                                const ClusteringParameters &parameters) noexcept -> bool {
        return cluster1.overlapsWithShortestSegmentFraction(
            cluster2, parameters.mergingStrandSpecificity, parameters.minSegmentOverlapFraction);
    };

    void greedyMerge(std::list<Interaction>::iterator seedIt) {
        bool additionalMerge = true;

        while (additionalMerge) {
            additionalMerge = false;

            for (auto iter = openClusterQueue.begin(); iter != openClusterQueue.end();) {
                if (iter == seedIt) {
                    ++iter;
                    continue;
                }

                if (clustersOverlap(*iter, *seedIt, parameters) &&
                    seedIt->merge(*iter, parameters.mergingStrandSpecificity)) {
                    iter = openClusterQueue.erase(iter);

                    additionalMerge = true;
                    break;
                }

                ++iter;
            }
        }
    }
};

}  // namespace pipelines::postprocess
