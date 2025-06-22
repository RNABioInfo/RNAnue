#pragma once

// Standard

// Internal
#include <algorithm>
#include <forward_list>
#include <functional>
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

            bool clusterMerged = false;
            auto prevIter = openClusterQueue.before_begin();

            // Try merging with clusters in the open queue
            for (auto iter = openClusterQueue.begin(); iter != openClusterQueue.end();) {
                if (clustersOverlap(*iter, cluster, parameters) &&
                    iter->merge(cluster, parameters.mergingStrandSpecificity)) {
                    clusterMerged = true;
                    break;
                }

                if (cluster.isBefore(*iter)) {
                    finishedClusters.emplace_back(std::move(*iter));
                    iter = openClusterQueue.erase_after(prevIter);
                    continue;
                }

                prevIter = iter;
                ++iter;
            }

            // If not merged, insert new cluster into the queue
            if (!clusterMerged) {
                openClusterQueue.emplace_after(prevIter, std::move(cluster));
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
    std::forward_list<Interaction> openClusterQueue;
    std::vector<Interaction> finishedClusters;

    static auto clustersOverlap(const Interaction &cluster1, const Interaction &cluster2,
                                const ClusteringParameters &parameters) noexcept -> bool {
        return cluster1.overlapsWithShortestSegmentFraction(
            cluster2, parameters.mergingStrandSpecificity, parameters.minSegmentOverlapFraction);
    };
};

}  // namespace pipelines::postprocess
