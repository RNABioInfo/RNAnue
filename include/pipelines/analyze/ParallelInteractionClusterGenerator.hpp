#pragma once

// Standard
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

// Internal
#include "ClusteringParameters.hpp"
#include "FeatureAnnotator.hpp"
#include "InteractionCluster.hpp"
#include "InteractionClusterGenerator.hpp"
#include "PartiallyAnnotatedInteractionCluster.hpp"

namespace pipelines::analyze {

using namespace annotation;

class ParallelInteractionClusterGenerator {
   public:
    ParallelInteractionClusterGenerator(std::shared_ptr<const FeatureAnnotator> featureAnnotator,
                                        ClusteringParameters parameters) noexcept
        : parameters(parameters), featureAnnotator(std::move(featureAnnotator)) {}

    struct Result {
        AnnotatedInteractionClusters annotatedClusters;
        FeatureCountsByFeatureID featureCounts;
        annotation::FeatureAnnotator supplementaryFeatureAnnotator;
    };

    /**
     * @brief Merges overlapping interaction clusters and returns these.
     *
     * This function takes a list of interaction clusters, sorts them from back to front,
     * and merges any overlapping clusters. This is done from back to front while closing clusters
     * that are further back than the current cluster.
     *
     * @param clusters A reference to the list of interaction clusters to be merged.
     * @return A list of finalized interaction clusters.
     */
    auto mergeClusters(std::vector<InteractionCluster>&& clusters, size_t threadCount,
                       size_t batchSize) -> Result;

   private:
    ClusteringParameters parameters;

    InteractionClusterGenerator::Result clusteringResults{};

    std::shared_ptr<const FeatureAnnotator> featureAnnotator;

    void annotatePartiallyAnnotatedClusters(
        const FeatureAnnotator& supplementaryFeatureAnnotator) noexcept;

    void logClusteringStatus() const noexcept;
};

}  // namespace pipelines::analyze
