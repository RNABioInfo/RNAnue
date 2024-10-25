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
#include "FeatureAnnotator.hpp"
#include "InteractionCluster.hpp"
#include "Orientation.hpp"
#include "PartiallyAnnotatedInteractionCluster.hpp"

namespace pipelines::analyze {

class InteractionClusterGenerator {
   public:
    InteractionClusterGenerator(std::string sampleName,
                                std::shared_ptr<annotation::FeatureAnnotator> featureAnnotator,
                                std::deque<std::string> referenceIDs,
                                annotation::Orientation featureOrientation, size_t minReadCount,
                                int graceDistance) noexcept
        : sampleName(std::move(sampleName)),
          featureAnnotator(std::move(featureAnnotator)),
          referenceIDs(std::move(referenceIDs)),
          featureOrientation(featureOrientation),
          minReadCount(minReadCount),
          graceDistance(graceDistance) {}

    using AnnotatedInteractionClusters = std::vector<AnnotatedInteractionCluster>;
    using FeatureCounts = std::unordered_map<std::string, size_t>;
    struct Result;

    auto mergeClusters(std::vector<InteractionCluster>&& clusters) -> Result;

   private:
    AnnotatedInteractionClusters finishedClusters;
    std::vector<PartiallyAnnotatedInteractionCluster> partiallyAnnotatedClusters;

    std::forward_list<InteractionCluster> openClusterQueue;

    std::string sampleName;

    std::shared_ptr<annotation::FeatureAnnotator> featureAnnotator;
    std::deque<std::string> referenceIDs;
    annotation::Orientation featureOrientation;

    annotation::FeatureAnnotator supplementaryFeatureAnnotator;

    size_t minReadCount;
    int graceDistance;

    FeatureCounts featureCounts;

    size_t includedClusterCount = 0;
    size_t excludedClusterCount = 0;

    void finalizeCluster(InteractionCluster&& cluster) noexcept;

    auto annotateCluster(InteractionCluster&& cluster) noexcept
        -> std::variant<AnnotatedInteractionCluster, PartiallyAnnotatedInteractionCluster>;

    void annotateSupplementaryFeatures() noexcept;

    static auto clusterIsBeforeOpenCluster(const InteractionCluster& cluster,
                                           const InteractionCluster& openCluster) noexcept -> bool;
    void logClusteringStatus() const noexcept;
};

struct InteractionClusterGenerator::Result {
    AnnotatedInteractionClusters annotatedClusters;
    FeatureCounts featureCounts;
    annotation::FeatureAnnotator supplementaryFeatureAnnotator;
};

}  // namespace pipelines::analyze
