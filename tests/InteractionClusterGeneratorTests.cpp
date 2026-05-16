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
#include "ClusteringParameters.hpp"
#include "FeatureAnnotator.hpp"
#include "GenomicFeature.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "InteractionCluster.hpp"
#include "InteractionClusterGenerator.hpp"
#include "Region.hpp"
#include "SortedGenomicRegionPair.hpp"

using namespace dataTypes;
using namespace pipelines::analyze;

namespace {

using ClusterShape = std::tuple<int, int, int, GenomicStrand, int, int, int, GenomicStrand, size_t>;

auto makeCluster(std::string recordID, int firstStart, int firstEnd, int secondStart,
                 int secondEnd, GenomicStrand firstStrand = GenomicStrand::FORWARD,
                 GenomicStrand secondStrand = GenomicStrand::FORWARD,
                 int firstReferenceID = 0, int secondReferenceID = 1) -> InteractionCluster {
    return InteractionCluster{
        SortedGenomicRegionPair{
            GenomicRegion{firstReferenceID, Region{.startPosition = firstStart,
                                                   .endPosition = firstEnd},
                          firstStrand},
            GenomicRegion{secondReferenceID,
                          Region{.startPosition = secondStart, .endPosition = secondEnd},
                          secondStrand}},
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

auto params(ClusteringMergeParameterVariant mergeParameter,
            GenomicStrandSpecificity strandSpecificity = GenomicStrandSpecificity::UNSPECIFIC)
    -> ClusteringParameters {
    return {.clusterMergeParameter = mergeParameter,
            .clusterMergingStrandSpecificity = strandSpecificity,
            .maxClusterSelfOverlapFraction = 1.0,
            .minimumClusterTrascriptContribution = 1.0F,
            .featureOrientation = GenomicOrientation::BOTH};
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

auto partialShapes(InteractionClusterGenerator::Result& result) -> std::vector<ClusterShape> {
    std::vector<ClusterShape> shapes;
    shapes.reserve(result.partiallyAnnotatedClusters.size());

    for (const auto& cluster : result.partiallyAnnotatedClusters) {
        shapes.emplace_back(shapeOf(cluster));
    }

    std::ranges::sort(shapes);
    return shapes;
}

}  // namespace

TEST(InteractionClusterGeneratorTests, BridgeReadConnectsPairwiseOverlappingComponents) {
    std::vector<InteractionCluster> clusters{
        makeCluster("left", 0, 10, 100, 110),
        makeCluster("bridge", 8, 18, 108, 118),
        makeCluster("right", 16, 26, 116, 126),
    };

    InteractionClusterGenerator generator{emptyFeatureAnnotator(),
                                          params(ClusterOverlapToleranceMergeParameter{0})};

    auto result = generator.mergeClusters(std::move(clusters));

    ASSERT_EQ(result.partiallyAnnotatedClusters.size(), 1UL);
    EXPECT_EQ(shapeOf(result.partiallyAnnotatedClusters.front()),
              ClusterShape(0, 0, 26, GenomicStrand::FORWARD, 1, 100, 126,
                           GenomicStrand::FORWARD, 3));
}

TEST(InteractionClusterGeneratorTests, RequiresBothChimericArmsToOverlap) {
    std::vector<InteractionCluster> clusters{
        makeCluster("first", 0, 10, 100, 110),
        makeCluster("second", 5, 15, 130, 140),
    };

    InteractionClusterGenerator generator{emptyFeatureAnnotator(),
                                          params(ClusterOverlapToleranceMergeParameter{0})};

    auto result = generator.mergeClusters(std::move(clusters));

    EXPECT_EQ(partialShapes(result),
              (std::vector<ClusterShape>{
                  {0, 0, 10, GenomicStrand::FORWARD, 1, 100, 110, GenomicStrand::FORWARD, 1},
                  {0, 5, 15, GenomicStrand::FORWARD, 1, 130, 140, GenomicStrand::FORWARD, 1},
              }));
}

TEST(InteractionClusterGeneratorTests, PositiveClusterDistanceMaterializesGappedComponent) {
    std::vector<InteractionCluster> clusters{
        makeCluster("left", 0, 10, 100, 110),
        makeCluster("right", 15, 25, 115, 125),
    };

    InteractionClusterGenerator generator{
        emptyFeatureAnnotator(),
        params(ClusterOverlapToleranceMergeParameter{clusterDistanceToOverlapTolerance(5)})};

    auto result = generator.mergeClusters(std::move(clusters));

    ASSERT_EQ(result.partiallyAnnotatedClusters.size(), 1UL);
    EXPECT_EQ(shapeOf(result.partiallyAnnotatedClusters.front()),
              ClusterShape(0, 0, 25, GenomicStrand::FORWARD, 1, 100, 125,
                           GenomicStrand::FORWARD, 2));
}

TEST(InteractionClusterGeneratorTests, LongOlderIntervalIsNotFinalizedBeforeAValidMerge) {
    std::vector<InteractionCluster> clusters{
        makeCluster("long", 10, 20, 100, 200),
        makeCluster("separate", 0, 5, 120, 150),
        makeCluster("overlaps-long", 15, 25, 90, 110),
    };

    InteractionClusterGenerator generator{emptyFeatureAnnotator(),
                                          params(ClusterOverlapToleranceMergeParameter{0})};

    auto result = generator.mergeClusters(std::move(clusters));

    EXPECT_EQ(partialShapes(result),
              (std::vector<ClusterShape>{
                  {0, 0, 5, GenomicStrand::FORWARD, 1, 120, 150, GenomicStrand::FORWARD, 1},
                  {0, 10, 25, GenomicStrand::FORWARD, 1, 90, 200, GenomicStrand::FORWARD, 2},
              }));
}

TEST(InteractionClusterGeneratorTests, ShortestSegmentFractionUsesBothArms) {
    std::vector<InteractionCluster> clusters{
        makeCluster("left", 0, 10, 100, 110),
        makeCluster("middle", 5, 15, 105, 115),
        makeCluster("right", 14, 24, 114, 124),
    };

    InteractionClusterGenerator generator{
        emptyFeatureAnnotator(), params(ShortestClusterOverlapFractionMergeParameter{0.5F})};

    auto result = generator.mergeClusters(std::move(clusters));

    EXPECT_EQ(partialShapes(result),
              (std::vector<ClusterShape>{
                  {0, 0, 15, GenomicStrand::FORWARD, 1, 100, 115, GenomicStrand::FORWARD, 2},
                  {0, 14, 24, GenomicStrand::FORWARD, 1, 114, 124, GenomicStrand::FORWARD, 1},
              }));
}

TEST(InteractionClusterGeneratorTests, StrandSpecificClusteringKeepsOppositeStrandsSeparate) {
    std::vector<InteractionCluster> clusters{
        makeCluster("forward", 0, 10, 100, 110, GenomicStrand::FORWARD),
        makeCluster("reverse", 5, 15, 105, 115, GenomicStrand::REVERSE),
    };

    InteractionClusterGenerator generator{
        emptyFeatureAnnotator(),
        params(ClusterOverlapToleranceMergeParameter{0}, GenomicStrandSpecificity::SPECIFIC)};

    auto result = generator.mergeClusters(std::move(clusters));

    EXPECT_EQ(partialShapes(result),
              (std::vector<ClusterShape>{
                  {0, 0, 10, GenomicStrand::FORWARD, 1, 100, 110, GenomicStrand::FORWARD, 1},
                  {0, 5, 15, GenomicStrand::REVERSE, 1, 105, 115, GenomicStrand::FORWARD, 1},
              }));
}

TEST(InteractionClusterGeneratorTests, UnspecificClusteringMergesAndCollapsesMixedStrands) {
    std::vector<InteractionCluster> clusters{
        makeCluster("forward", 0, 10, 100, 110, GenomicStrand::FORWARD),
        makeCluster("reverse", 5, 15, 105, 115, GenomicStrand::REVERSE),
    };

    InteractionClusterGenerator generator{
        emptyFeatureAnnotator(),
        params(ClusterOverlapToleranceMergeParameter{0}, GenomicStrandSpecificity::UNSPECIFIC)};

    auto result = generator.mergeClusters(std::move(clusters));

    ASSERT_EQ(result.partiallyAnnotatedClusters.size(), 1UL);
    EXPECT_EQ(shapeOf(result.partiallyAnnotatedClusters.front()),
              ClusterShape(0, 0, 15, GenomicStrand::NONE, 1, 100, 115,
                           GenomicStrand::FORWARD, 2));
}

TEST(InteractionClusterGeneratorTests, DifferentReferencePairsCannotMerge) {
    std::vector<InteractionCluster> clusters{
        makeCluster("first-reference-pair", 0, 10, 100, 110, GenomicStrand::FORWARD,
                    GenomicStrand::FORWARD, 0, 1),
        makeCluster("second-reference-pair", 5, 15, 105, 115, GenomicStrand::FORWARD,
                    GenomicStrand::FORWARD, 0, 2),
    };

    InteractionClusterGenerator generator{emptyFeatureAnnotator(),
                                          params(ClusterOverlapToleranceMergeParameter{0})};

    auto result = generator.mergeClusters(std::move(clusters));

    EXPECT_EQ(partialShapes(result),
              (std::vector<ClusterShape>{
                  {0, 0, 10, GenomicStrand::FORWARD, 1, 100, 110, GenomicStrand::FORWARD, 1},
                  {0, 5, 15, GenomicStrand::FORWARD, 2, 105, 115, GenomicStrand::FORWARD, 1},
              }));
}

// NOLINTEND(readability-magic-numbers)
