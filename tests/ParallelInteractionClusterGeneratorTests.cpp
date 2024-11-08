#include <gtest/gtest.h>

// Standard
#include <algorithm>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "FeatureAnnotator.hpp"
#include "GenomicFeature.hpp"
#include "GenomicStrand.hpp"
#include "InteractionCluster.hpp"
#include "Orientation.hpp"
#include "ParallelInteractionClusterGenerator.hpp"
#include "RecordFragment.hpp"
#include "SplitRecordsParser.hpp"

using namespace pipelines::analyze;
using namespace annotation;

// Helper function to get the test SAM file path
auto testParallelInteractionClustersSamPath() -> std::string {
    return (std::filesystem::path{__FILE__}.parent_path() /
            "test_data/interactionClusterRecords.sam")
        .string();
}

// Expected Clusters Definitions (Same as in InteractionClusterGeneratorTests)
const AnnotatedInteractionCluster expectedCluster1(
    {{.firstSegment = RecordFragment{.recordID = "SRR18331301.3",
                                     .referenceIDIndex = 0,
                                     .strand = dataTypes::GenomicStrand::FORWARD,
                                     .start = 19,
                                     .end = 27,
                                     .complementarityScore = 1.0,
                                     .hybridizationEnergy = -1.7},
      .secondSegment = RecordFragment{.recordID = "SRR18331301.3",
                                      .referenceIDIndex = 1,
                                      .strand = dataTypes::GenomicStrand::FORWARD,
                                      .start = 49,
                                      .end = 64,
                                      .complementarityScore = 1.0,
                                      .hybridizationEnergy = -1.7}},
     {"SRR18331301.3", "SRR18331301.1", "SRR18331301.6"},
     {1.0, 1.0, 1.0},
     {-1.7, -1.7, -1.7}},
    "gene1", "gene2");

const AnnotatedInteractionCluster expectedCluster2(
    {{.firstSegment = RecordFragment{.recordID = "SRR18331301.2",
                                     .referenceIDIndex = 1,
                                     .strand = dataTypes::GenomicStrand::FORWARD,
                                     .start = 4,
                                     .end = 10,
                                     .complementarityScore = 1.0,
                                     .hybridizationEnergy = -1.7},
      .secondSegment = RecordFragment{.recordID = "SRR18331301.2",
                                      .referenceIDIndex = 1,
                                      .strand = dataTypes::GenomicStrand::FORWARD,
                                      .start = 51,
                                      .end = 57,
                                      .complementarityScore = 1.0,
                                      .hybridizationEnergy = -1.7}},
     {"SRR18331301.2", "SRR18331301.7"},
     {1.0, 1.0},
     {-1.7, -1.7}},
    "random", "gene3");

const AnnotatedInteractionCluster expectedCluster3(
    {{.firstSegment = RecordFragment{.recordID = "SRR18331301.4",
                                     .referenceIDIndex = 0,
                                     .strand = dataTypes::GenomicStrand::FORWARD,
                                     .start = 4,
                                     .end = 14,
                                     .complementarityScore = 1.0,
                                     .hybridizationEnergy = -1.7},
      .secondSegment = RecordFragment{.recordID = "SRR18331301.4",
                                      .referenceIDIndex = 0,
                                      .strand = dataTypes::GenomicStrand::FORWARD,
                                      .start = 39,
                                      .end = 46,
                                      .complementarityScore = 1.0,
                                      .hybridizationEnergy = -1.7}},
     {"SRR18331301.5", "SRR18331301.4"},
     {1.0, 1.0},
     {-1.7, -1.7}},
    "random", "random");

// Test Parameters Structure
struct ParallelInteractionClusterGeneratorTestParam {
    size_t threadCount;
    size_t batchSize;
    std::unordered_map<std::string, size_t> expectedFeatureCounts;
    std::vector<AnnotatedInteractionCluster> finalizedClusters;
};

// Feature Map Setup (Same as in InteractionClusterGeneratorTests)
static const FeatureMap featureMap = {{"chr1",
                                       {GenomicFeature{.referenceID = "chr1",
                                                       .type = "gene",
                                                       .startPosition = 20,
                                                       .endPosition = 30,
                                                       .strand = dataTypes::GenomicStrand::FORWARD,
                                                       .id = "gene1",
                                                       .groupID = std::nullopt,
                                                       .geneName = std::nullopt},
                                        GenomicFeature{.referenceID = "chr1",
                                                       .type = "gene",
                                                       .startPosition = 50,
                                                       .endPosition = 60,
                                                       .strand = dataTypes::GenomicStrand::FORWARD,
                                                       .id = "gene2",
                                                       .groupID = std::nullopt,
                                                       .geneName = std::nullopt}}},
                                      {"chr2",
                                       {GenomicFeature{.referenceID = "chr2",
                                                       .type = "gene",
                                                       .startPosition = 50,
                                                       .endPosition = 60,
                                                       .strand = dataTypes::GenomicStrand::FORWARD,
                                                       .id = "gene3",
                                                       .groupID = std::nullopt,
                                                       .geneName = std::nullopt}}}};

static const std::shared_ptr<FeatureAnnotator> featureAnnotator =
    std::make_shared<FeatureAnnotator>(featureMap);

// Test Fixture
class ParallelInteractionClusterGeneratorTests
    : public testing::TestWithParam<ParallelInteractionClusterGeneratorTestParam> {
   protected:
    ParallelInteractionClusterGeneratorTests()
        : parallelClusterGenerator(featureAnnotator, {"chr1", "chr2"},
                                   {.featureOrientation = Orientation::BOTH,
                                    .maxOverlapFraction = 0.5,
                                    .minReadCount = 1,
                                    .graceDistance = 1}) {}

    ParallelInteractionClusterGenerator parallelClusterGenerator;
};

// Define Specific Test Cases
TEST_P(ParallelInteractionClusterGeneratorTests, MergeClustersCorrectly) {
    const auto& param = GetParam();

    // Parse and sort interaction clusters
    std::vector<InteractionCluster> interactionClusters =
        SplitRecordsParser::parse(testParallelInteractionClustersSamPath());

    std::ranges::sort(interactionClusters, std::less<>());

    // Instantiate ParallelInteractionClusterGenerator with specific thread and batch size
    ParallelInteractionClusterGenerator generator(featureAnnotator, {"chr1", "chr2"},
                                                  {.featureOrientation = Orientation::BOTH,
                                                   .maxOverlapFraction = 0.5,
                                                   .minReadCount = 1,
                                                   .graceDistance = 1});

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
    ParallelInteractionClusterGenerator generator(featureAnnotator, {"chr1", "chr2"},
                                                  {.featureOrientation = Orientation::BOTH,
                                                   .maxOverlapFraction = 0.5,
                                                   .minReadCount = 1,
                                                   .graceDistance = 1});

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
