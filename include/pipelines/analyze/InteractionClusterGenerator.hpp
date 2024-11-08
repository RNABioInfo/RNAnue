#pragma once

// Standard
#include <cstddef>
#include <deque>
#include <forward_list>
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
#include "GenomicRegion.hpp"
#include "InteractionCluster.hpp"
#include "PartiallyAnnotatedInteractionCluster.hpp"

namespace pipelines::analyze {

class InteractionClusterGenerator {
   public:
    using AnnotatedInteractionClusters = std::vector<AnnotatedInteractionCluster>;
    using FeatureCounts = std::unordered_map<std::string, size_t>;

    struct Result {
        AnnotatedInteractionClusters finishedClusters;
        std::vector<PartiallyAnnotatedInteractionCluster> partiallyAnnotatedClusters;
        std::vector<GenomicRegion> supplementaryFeatureRegions;
        FeatureCounts featureCounts;
        size_t includedClusterCount;
        size_t excludedClusterCount;
    };

    InteractionClusterGenerator(std::shared_ptr<annotation::FeatureAnnotator> featureAnnotator,
                                std::deque<std::string> referenceIDs,
                                ClusteringParameters parameters) noexcept
        : parameters(parameters),
          featureAnnotator(std::move(featureAnnotator)),
          referenceIDs(std::move(referenceIDs)) {}

    auto mergeClusters(std::vector<InteractionCluster>&& clusters) -> Result;

   private:
    ClusteringParameters parameters;

    AnnotatedInteractionClusters finishedClusters;
    std::vector<PartiallyAnnotatedInteractionCluster> partiallyAnnotatedClusters;

    std::vector<GenomicRegion> supplementaryFeatureRegions;

    std::forward_list<InteractionCluster> openClusterQueue;

    std::shared_ptr<annotation::FeatureAnnotator> featureAnnotator;
    std::deque<std::string> referenceIDs;

    FeatureCounts featureCounts;

    size_t includedClusterCount = 0;
    size_t excludedClusterCount = 0;

    void finalizeCluster(InteractionCluster&& cluster) noexcept;

    auto annotateCluster(InteractionCluster&& cluster) noexcept
        -> std::variant<AnnotatedInteractionCluster, PartiallyAnnotatedInteractionCluster>;
};

}  // namespace pipelines::analyze
