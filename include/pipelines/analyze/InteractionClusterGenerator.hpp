#pragma once

// Standard
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "ClusteringParameters.hpp"
#include "FeatureAnnotator.hpp"
#include "GenomicFeature.hpp"
#include "InteractionCluster.hpp"
#include "PartiallyAnnotatedInteractionCluster.hpp"

namespace pipelines::analyze {

using AnnotatedInteractionClusters = std::vector<AnnotatedInteractionCluster>;
using PartiallyAnnotatedClusters = std::vector<PartiallyAnnotatedInteractionCluster>;
using FeatureCountsByFeatureID = std::unordered_map<std::string, float>;

using namespace annotation;

class InteractionClusterGenerator {
   public:
    InteractionClusterGenerator(std::shared_ptr<const FeatureAnnotator> featureAnnotator,
                                ClusteringParameters parameters) noexcept
        : parameters(parameters), featureAnnotator(std::move(featureAnnotator)) {}

    struct Result {
        AnnotatedInteractionClusters finishedClusters;
        PartiallyAnnotatedClusters partiallyAnnotatedClusters;
        FeatureMap supplementaryFeatureMap;
        FeatureCountsByFeatureID featureCounts;
        size_t includedClusterCount{0};
        size_t excludedClusterCount{0};

        void merge(Result&& other) noexcept;

        [[nodiscard]] constexpr auto totalClusterCount() const noexcept -> size_t {
            return includedClusterCount + excludedClusterCount;
        };
    };

    /**
     * @brief Merges overlapping interaction clusters and returns these.
     *
     * This function takes a sorted list of interaction clusters,
     * and merges any overlapping clusters. This is done from back to front while closing clusters
     * that are further back than the current cluster.
     *
     * @param clusters A forwarding reference to the sorted list of interaction clusters to be
     * merged.
     * @return A list of finalized interaction clusters.
     */
    auto mergeClusters(std::vector<InteractionCluster>&& clusters) -> Result;

    auto mergeClusterGroup(std::vector<InteractionCluster>&& clusters) -> Result;

   private:
    ClusteringParameters parameters;

    AnnotatedInteractionClusters finishedClusters;
    std::vector<PartiallyAnnotatedInteractionCluster> partiallyAnnotatedClusters;

    FeatureMap supplementaryFeatureRegions;

    std::shared_ptr<const FeatureAnnotator> featureAnnotator;

    FeatureCountsByFeatureID featureCountsByFeatureID;

    size_t includedClusterCount = 0;
    size_t excludedClusterCount = 0;

    [[nodiscard]] static auto clustersOverlap(const InteractionCluster& cluster1,
                                              const InteractionCluster& cluster2,
                                              const ClusteringParameters& parameters) noexcept
        -> bool;

    [[nodiscard]] auto clusterPassesFilters(const InteractionCluster& cluster) const noexcept
        -> bool;

    [[nodiscard]] auto annotateCluster(InteractionCluster&& cluster) noexcept
        -> std::variant<AnnotatedInteractionCluster, PartiallyAnnotatedInteractionCluster>;

    void attributeCluster(AnnotatedInteractionCluster&& cluster) noexcept;

    void attributeCluster(PartiallyAnnotatedInteractionCluster&& cluster) noexcept;

    void finalizeCluster(InteractionCluster&& cluster) noexcept;

    void finalizeMergedClusters(std::vector<InteractionCluster>&& clusters) noexcept;

    [[nodiscard]] auto releaseResult() noexcept -> Result;
};

}  // namespace pipelines::analyze
