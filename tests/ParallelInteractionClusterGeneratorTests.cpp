// NOLINTBEGIN(readability-magic-numbers)

#include <gtest/gtest.h>

// Standard
#include <algorithm>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "ClusteringParameters.hpp"
#include "FeatureAnnotator.hpp"
#include "GenomicFeature.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "InteractionCluster.hpp"
#include "InteractionClusterGenerator.hpp"
#include "ParallelInteractionClusterGenerator.hpp"
#include "Region.hpp"
#include "SortedGenomicRegionPair.hpp"

using namespace dataTypes;
using namespace pipelines::analyze;

namespace {

using ClusterShape = std::tuple<int, int, int, GenomicStrand, int, int, int, GenomicStrand, size_t>;

auto makeCluster(std::string recordID, int firstStart, int firstEnd, int secondStart,
                 int secondEnd, int firstReferenceID = 0, int secondReferenceID = 1)
    -> InteractionCluster {
    return InteractionCluster{
        SortedGenomicRegionPair{
            GenomicRegion{firstReferenceID, Region{.startPosition = firstStart,
                                                   .endPosition = firstEnd},
                          GenomicStrand::FORWARD},
            GenomicRegion{secondReferenceID,
                          Region{.startPosition = secondStart, .endPosition = secondEnd},
                          GenomicStrand::FORWARD}},
        std::move(recordID),
        1.0,
        -1.0,
        1,
        1.0F};
}

auto emptyFeatureAnnotator() -> std::shared_ptr<const annotation::FeatureAnnotator> {
    const FeatureMap featureMap = {{0, {}}, {1, {}}, {2, {}}};
    return std::make_shared<const annotation::FeatureAnnotator>(featureMap);
}

auto params() -> ClusteringParameters {
    return {.clusterMergeParameter = ClusterOverlapToleranceMergeParameter{0},
            .clusterMergingStrandSpecificity = GenomicStrandSpecificity::UNSPECIFIC,
            .maxClusterSelfOverlapFraction = 1.0,
            .minimumClusterTrascriptContribution = 1.0F,
            .featureOrientation = GenomicOrientation::BOTH};
}

auto inputClusters() -> std::vector<InteractionCluster> {
    return {
        makeCluster("component-a-1", 0, 10, 100, 110),
        makeCluster("component-b-1", 200, 210, 500, 510),
        makeCluster("component-a-2", 8, 18, 108, 118),
        makeCluster("component-b-2", 208, 218, 508, 518),
        makeCluster("other-reference-pair", 8, 18, 108, 118, 0, 2),
    };
}

auto shapeOf(const InteractionCluster& cluster) -> ClusterShape {
    return {cluster.getFirstSegment().getReferenceIDIndex(),
            cluster.getFirstSegment().getStart(),
            cluster.getFirstSegment().getEnd(),
            cluster.getFirstSegment().getStrand(),
            cluster.getSecondSegment().getReferenceIDIndex(),
            cluster.getSecondSegment().getStart(),
            cluster.getSecondSegment().getEnd(),
            cluster.getSecondSegment().getStrand(),
            cluster.fragmentCount()};
}

auto shapesFrom(const std::vector<AnnotatedInteractionCluster>& clusters)
    -> std::vector<ClusterShape> {
    std::vector<ClusterShape> shapes;
    shapes.reserve(clusters.size());

    for (const auto& cluster : clusters) {
        shapes.emplace_back(shapeOf(cluster));
    }

    std::ranges::sort(shapes);
    return shapes;
}

auto expectedShapes() -> std::vector<ClusterShape> {
    std::vector<ClusterShape> shapes{
        {0, 0, 18, GenomicStrand::FORWARD, 1, 100, 118, GenomicStrand::FORWARD, 2},
        {0, 8, 18, GenomicStrand::FORWARD, 2, 108, 118, GenomicStrand::FORWARD, 1},
        {0, 200, 218, GenomicStrand::FORWARD, 1, 500, 518, GenomicStrand::FORWARD, 2},
    };
    std::ranges::sort(shapes);
    return shapes;
}

}  // namespace

TEST(ParallelInteractionClusterGeneratorTests, SingleThreadMatchesConnectedComponentResult) {
    ParallelInteractionClusterGenerator generator{emptyFeatureAnnotator(), params()};

    auto result = generator.mergeClusters(inputClusters(), 1, 1);

    EXPECT_EQ(shapesFrom(result.annotatedClusters), expectedShapes());
}

TEST(ParallelInteractionClusterGeneratorTests, MultipleThreadsMatchSingleThreadRegardlessOfBatchSize) {
    ParallelInteractionClusterGenerator singleThreadGenerator{emptyFeatureAnnotator(), params()};
    ParallelInteractionClusterGenerator parallelGenerator{emptyFeatureAnnotator(), params()};

    auto singleThreadResult = singleThreadGenerator.mergeClusters(inputClusters(), 1, 1);
    auto parallelResult = parallelGenerator.mergeClusters(inputClusters(), 4, 1);

    EXPECT_EQ(shapesFrom(parallelResult.annotatedClusters),
              shapesFrom(singleThreadResult.annotatedClusters));
}

TEST(ParallelInteractionClusterGeneratorTests, HandlesEmptyInput) {
    ParallelInteractionClusterGenerator generator{emptyFeatureAnnotator(), params()};

    auto result = generator.mergeClusters({}, 4, 1);

    EXPECT_TRUE(result.annotatedClusters.empty());
    EXPECT_TRUE(result.featureCounts.empty());
}

// NOLINTEND(readability-magic-numbers)
