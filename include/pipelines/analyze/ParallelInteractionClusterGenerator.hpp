#pragma once

// Standard
#include <string>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "ClusteringParameters.hpp"
#include "FeatureAnnotator.hpp"
#include "InteractionClusterGenerator.hpp"

namespace pipelines::analyze {

class ParallelInteractionClusterGenerator {
   public:
    using AnnotatedInteractionClusters = std::vector<AnnotatedInteractionCluster>;
    using FeatureCounts = std::unordered_map<std::string, size_t>;

    struct Result {
        AnnotatedInteractionClusters annotatedClusters;
        FeatureCounts featureCounts;
        annotation::FeatureAnnotator supplementaryFeatureAnnotator;
    };

    ParallelInteractionClusterGenerator(
        std::shared_ptr<annotation::FeatureAnnotator> featureAnnotator,
        std::deque<std::string> referenceIDs, ClusteringParameters parameters) noexcept
        : parameters(parameters),
          featureAnnotator(std::move(featureAnnotator)),
          referenceIDs(std::move(referenceIDs)) {}

    auto mergeClusters(std::vector<InteractionCluster>&& clusters, size_t threadCount,
                       size_t batchSize) -> Result;

   private:
    ClusteringParameters parameters;

    AnnotatedInteractionClusters finishedClusters;
    std::vector<PartiallyAnnotatedInteractionCluster> partiallyAnnotatedClusters;

    std::shared_ptr<annotation::FeatureAnnotator> featureAnnotator;
    std::deque<std::string> referenceIDs;

    annotation::FeatureAnnotator supplementaryFeatureAnnotator;

    FeatureCounts featureCounts;

    size_t includedClusterCount = 0;
    size_t excludedClusterCount = 0;

    void annotateSupplementaryFeatures() noexcept;

    void mergeClusteringResults(InteractionClusterGenerator::Result&& result) noexcept;

    void logClusteringStatus() const noexcept;
};

}  // namespace pipelines::analyze
