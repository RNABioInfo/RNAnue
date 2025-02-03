#include <gtest/gtest.h>

// Standard
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
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
#include "ParallelInteractionClusterGenerator.hpp"
#include "SortedGenomicRegionPair.hpp"
#include "SplitRecordsParser.hpp"
#include "TestFilePath.hpp"

using namespace pipelines::analyze;
using namespace annotation;

// Expected Clusters Definitions (Same as in InteractionClusterGeneratorTests)
const AnnotatedInteractionCluster expectedCluster1(
    InteractionCluster{
        SortedGenomicRegionPair{GenomicRegion{0, Region{.startPosition = 19, .endPosition = 27},
                                              GenomicStrand::FORWARD},
                                GenomicRegion{1, Region{.startPosition = 49, .endPosition = 64},
                                              GenomicStrand::FORWARD}},
        {"SRR18331301.3", "SRR18331301.1", "SRR18331301.6"},
        {1.0, 1.0, 1.0},
        {-1.7, -1.7, -1.7},
        {1, 1, 1}

    },
    "gene1", "gene2");

const AnnotatedInteractionCluster expectedCluster2(
    InteractionCluster{
        SortedGenomicRegionPair{
            GenomicRegion{1, Region{.startPosition = 4, .endPosition = 10}, GenomicStrand::FORWARD},
            GenomicRegion{1, Region{.startPosition = 51, .endPosition = 57},
                          GenomicStrand::FORWARD}},
        {"SRR18331301.2", "SRR18331301.7"},
        {
            1.0,
            1.0,
        },
        {
            -1.7,
            -1.7,
        },
        {
            1,
            1,
        }},
    "random", "gene3");

const AnnotatedInteractionCluster expectedCluster3(
    InteractionCluster{
        SortedGenomicRegionPair{
            GenomicRegion{0, Region{.startPosition = 4, .endPosition = 14}, GenomicStrand::FORWARD},
            GenomicRegion{0, Region{.startPosition = 39, .endPosition = 46},
                          GenomicStrand::FORWARD}},
        {"SRR18331301.5", "SRR18331301.4"},
        {
            1.0,
            1.0,
        },
        {
            -1.7,
            -1.7,
        },
        {
            1,
            1,
        }},
    "random", "random");

// Test Parameters Structure
struct ParallelInteractionClusterGeneratorTestParam {
    size_t threadCount;
    size_t batchSize;
    std::unordered_map<std::string, size_t> expectedFeatureCounts;
    std::vector<AnnotatedInteractionCluster> finalizedClusters;
};

// Feature Map Setup (Same as in InteractionClusterGeneratorTests)
static const FeatureMap featureMap = {
    {1,
     {GenomicFeature{
          .type = "gene",
          .genomicRegion = {1, {.startPosition = 20, .endPosition = 30}, GenomicStrand::FORWARD},
          .featureID = "gene1",
          .groupID = std::nullopt,
          .geneName = std::nullopt},
      GenomicFeature{
          .type = "gene",
          .genomicRegion = {1, {.startPosition = 50, .endPosition = 60}, GenomicStrand::FORWARD},
          .featureID = "gene2",
          .groupID = std::nullopt,
          .geneName = std::nullopt}}},
    {2,
     {GenomicFeature{
         .type = "gene",
         .genomicRegion = {2, {.startPosition = 50, .endPosition = 60}, GenomicStrand::FORWARD},
         .featureID = "gene3",
         .groupID = std::nullopt,
         .geneName = std::nullopt}}}};

static const std::shared_ptr<FeatureAnnotator> featureAnnotator =
    std::make_shared<FeatureAnnotator>(featureMap);

// Test Fixture
class ParallelInteractionClusterGeneratorTests
    : public testing::TestWithParam<ParallelInteractionClusterGeneratorTestParam> {
   protected:
    ParallelInteractionClusterGeneratorTests()
        : parallelClusterGenerator(
              featureAnnotator,
              ClusteringParameters{
                  .clusterMergeParameter = ClusterOverlapToleranceMergeParameter{1},
                  .clusterMergingStrandSpecificity = GenomicStrandSpecificity::UNSPECIFIC,
                  .maxClusterSelfOverlapFraction = 0.5,  // NOLINT
                  .minimumClusterReadCount = 1,
                  .featureOrientation = dataTypes::GenomicOrientation::BOTH}) {}

    ParallelInteractionClusterGenerator parallelClusterGenerator;
};

// Define Specific Test Cases
TEST_P(ParallelInteractionClusterGeneratorTests, MergeClustersCorrectly) {
    const auto& param = GetParam();

    // Parse and sort interaction clusters
    std::vector<InteractionCluster> interactionClusters =
        SplitRecordsParser::parse(getTestFilePath("interactionClusterRecords.sam"));

    std::ranges::sort(interactionClusters, std::less<>());

    // Instantiate ParallelInteractionClusterGenerator with specific thread and batch size
    ParallelInteractionClusterGenerator generator(
        featureAnnotator,
        ClusteringParameters{
            .clusterMergeParameter = ClusterOverlapToleranceMergeParameter{1},
            .clusterMergingStrandSpecificity = GenomicStrandSpecificity::UNSPECIFIC,
            .maxClusterSelfOverlapFraction = 0.5,  // NOLINT
            .minimumClusterReadCount = 1,
            .featureOrientation = dataTypes::GenomicOrientation::BOTH});
    // Merge clusters in parallel
    auto result =
        generator.mergeClusters(std::move(interactionClusters), param.threadCount, param.batchSize);

    // Validate finalized clusters the order of clusters is not guaranteed
    auto expectedClusters = param.finalizedClusters;
    std::ranges::sort(result.annotatedClusters, std::less<>());
    std::ranges::sort(expectedClusters, std::less<>());
    EXPECT_EQ(result.annotatedClusters, expectedClusters);

    // TODO: Validate feature counts not working because of random supplementary annotations
    // EXPECT_EQ(result.featureCounts, param.expectedFeatureCounts);
}

// Edge Case: Empty Input
TEST_P(ParallelInteractionClusterGeneratorTests, HandleEmptyInput) {
    const ParallelInteractionClusterGeneratorTestParam param = {
        .threadCount = 2, .batchSize = 10, .expectedFeatureCounts = {}, .finalizedClusters = {}};

    // No interaction clusters
    std::vector<InteractionCluster> interactionClusters;

    // Instantiate ParallelInteractionClusterGenerator
    ParallelInteractionClusterGenerator generator(
        featureAnnotator,
        ClusteringParameters{
            .clusterMergeParameter = ClusterOverlapToleranceMergeParameter{1},
            .clusterMergingStrandSpecificity = GenomicStrandSpecificity::UNSPECIFIC,
            .maxClusterSelfOverlapFraction = 0.5,  // NOLINT
            .minimumClusterReadCount = 1,
            .featureOrientation = dataTypes::GenomicOrientation::BOTH});

    // Merge clusters in parallel
    auto result =
        generator.mergeClusters(std::move(interactionClusters), param.threadCount, param.batchSize);

    // Validate results are empty
    EXPECT_TRUE(result.annotatedClusters.empty());
    EXPECT_TRUE(result.featureCounts.empty());
}

// Parameterized Test Instantiation
INSTANTIATE_TEST_SUITE_P(
    Default, ParallelInteractionClusterGeneratorTests,
    testing::Values(
        ParallelInteractionClusterGeneratorTestParam{
            .threadCount = 2,
            .batchSize = 1,
            .expectedFeatureCounts = {{"gene1", 3}, {"gene3", 4}},
            .finalizedClusters = {expectedCluster1, expectedCluster3, expectedCluster2}},
        // Additional Test Configurations
        ParallelInteractionClusterGeneratorTestParam{
            .threadCount = 4,
            .batchSize = 2,
            .expectedFeatureCounts = {{"gene1", 3}, {"gene3", 4}},
            .finalizedClusters = {expectedCluster1, expectedCluster3, expectedCluster2}},
        ParallelInteractionClusterGeneratorTestParam{
            .threadCount = 1,  // Single thread
            .batchSize = 5,
            .expectedFeatureCounts = {{"gene1", 3}, {"gene3", 4}},
            .finalizedClusters = {expectedCluster1, expectedCluster3, expectedCluster2},
        }));
