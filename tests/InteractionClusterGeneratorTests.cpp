#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "AnnotatedInteractionCluster.hpp"
#include "FeatureAnnotator.hpp"
#include "GenomicFeature.hpp"
#include "GenomicStrand.hpp"
#include "InteractionCluster.hpp"
#include "InteractionClusterGenerator.hpp"
#include "Orientation.hpp"
#include "PartiallyAnnotatedInteractionCluster.hpp"
#include "RecordFragment.hpp"
#include "SplitRecordsParser.hpp"
#include "TestFilePath.hpp"

using namespace pipelines::analyze;

// Note: This is the expected clustering result:
// Cluster 1: 1,3,6 // Fully annotated; First segment: InteractionSegment(referenceIDIndex: 0,
// strand: +, start: 19, end: 27), Second segment id: InteractionSegment(referenceIDIndex: 1,
// strand: +, start: 49, end: 64)
//
// Cluster 2: 2,7 // Partially annotated; First segment: InteractionSegment(referenceIDIndex: 1,
// strand: +, start: 4, end: 10), Second segment id: InteractionSegment(referenceIDIndex: 1, strand:
// +, start: 51, end: 57)
//
// Cluster 3: 4,5 // Non annotated; First segment: InteractionSegment(referenceIDIndex: 0, strand:
// +, start: 4, end: 14), Second segment id: InteractionSegment(referenceIDIndex: 0, strand: +,
// start: 39, end: 46)

using namespace pipelines::analyze;

struct InteractionClusterGeneratorTestParam {
    std::unordered_map<std::string, size_t> expectedFeatureCounts;
    std::vector<AnnotatedInteractionCluster> finalizedClusters;
    std::vector<PartiallyAnnotatedInteractionCluster> partiallyAnnotatedClusters;
};

static const FeatureMap featureMap = {{"chr1",
                                       {GenomicFeature{.referenceID = "chr1",
                                                       .type = "gene",
                                                       .startPosition = 20,
                                                       .endPosition = 30,
                                                       .strand = GenomicStrand::FORWARD,
                                                       .id = "gene1",
                                                       .groupID = std::nullopt,
                                                       .geneName = std::nullopt},
                                        GenomicFeature{.referenceID = "chr1",
                                                       .type = "gene",
                                                       .startPosition = 50,
                                                       .endPosition = 60,
                                                       .strand = GenomicStrand::FORWARD,
                                                       .id = "gene2",
                                                       .groupID = std::nullopt,
                                                       .geneName = std::nullopt}}},
                                      {"chr2",
                                       {GenomicFeature{.referenceID = "chr2",
                                                       .type = "gene",
                                                       .startPosition = 50,
                                                       .endPosition = 60,
                                                       .strand = GenomicStrand::FORWARD,
                                                       .id = "gene3",
                                                       .groupID = std::nullopt,
                                                       .geneName = std::nullopt}}}};

static const std::shared_ptr<annotation::FeatureAnnotator> featureAnnotator =
    std::make_shared<annotation::FeatureAnnotator>(featureMap);

class InteractionClusterGeneratorTests
    : public testing::TestWithParam<InteractionClusterGeneratorTestParam> {
   protected:
    InteractionClusterGeneratorTests()
        : interactionClusterGenerator(featureAnnotator, {"chr1", "chr2"},
                                      {.featureOrientation = annotation::Orientation::BOTH,
                                       .maxOverlapFraction = 0.5,
                                       .minReadCount = 1,
                                       .graceDistance = 1}) {}

    InteractionClusterGenerator interactionClusterGenerator;
};

// TODO: Fix this test and use == operator for InteractionCluster2
TEST_P(InteractionClusterGeneratorTests, SplitRecordsAreSortedCorrectly) {
    const auto& param = GetParam();

    std::vector<InteractionCluster> interactionClusters =
        SplitRecordsParser::parse(getTestFilePath("interactionClusterRecords.sam"));

    std::ranges::sort(interactionClusters, std::less<>());

    EXPECT_EQ(interactionClusters.size(), 7U);

    auto result = interactionClusterGenerator.mergeClusters(std::move(interactionClusters));

    EXPECT_EQ(result.finishedClusters, param.finalizedClusters);
    EXPECT_EQ(result.partiallyAnnotatedClusters, param.partiallyAnnotatedClusters);
    EXPECT_EQ(result.featureCounts, param.expectedFeatureCounts);
}

const AnnotatedInteractionCluster cluster1(
    {{.firstSegment = RecordFragment{.recordID = "SRR18331301.3",
                                     .referenceIDIndex = 0,
                                     .strand = dataTypes::GenomicStrand::FORWARD,
                                     .start = 19,
                                     .end = 27,
                                     .complementarityScore = 1,
                                     .hybridizationEnergy = -1.7,
                                     .crosslinkingSiteCount = 1},
      .secondSegment = RecordFragment{.recordID = "SRR18331301.3",
                                      .referenceIDIndex = 1,
                                      .strand = dataTypes::GenomicStrand::FORWARD,
                                      .start = 49,
                                      .end = 64,
                                      .complementarityScore = 1,
                                      .hybridizationEnergy = -1.7,
                                      .crosslinkingSiteCount = 1}},
     {"SRR18331301.3", "SRR18331301.1", "SRR18331301.6"},
     {1.0, 1.0, 1.0},
     {-1.7, -1.7, -1.7},
     {1, 1, 1}},
    "gene1", "gene2");

const PartiallyAnnotatedInteractionCluster cluster2(
    {InteractionSegmentPair{
         .firstSegment = RecordFragment{.recordID = "SRR18331301.2",
                                        .referenceIDIndex = 1,
                                        .strand = dataTypes::GenomicStrand::FORWARD,
                                        .start = 4,
                                        .end = 10,
                                        .complementarityScore = 1.0,
                                        .hybridizationEnergy = -1.7,
                                        .crosslinkingSiteCount = 1},
         .secondSegment = RecordFragment{.recordID = "SRR18331301.2",
                                         .referenceIDIndex = 1,
                                         .strand = dataTypes::GenomicStrand::FORWARD,
                                         .start = 51,
                                         .end = 57,
                                         .complementarityScore = 1.0,
                                         .hybridizationEnergy = -1.7,
                                         .crosslinkingSiteCount = 1}},
     {"SRR18331301.2", "SRR18331301.7"},
     {1.0, 1.0},
     {-1.7, -1.7},
     {1, 1}},
    std::nullopt, "gene3");

const PartiallyAnnotatedInteractionCluster cluster3(
    {{.firstSegment = RecordFragment{.recordID = "SRR18331301.4",
                                     .referenceIDIndex = 0,
                                     .strand = dataTypes::GenomicStrand::FORWARD,
                                     .start = 4,
                                     .end = 14,
                                     .complementarityScore = 1.0,
                                     .hybridizationEnergy = -1.7,
                                     .crosslinkingSiteCount = 1},
      .secondSegment = RecordFragment{.recordID = "SRR18331301.4",
                                      .referenceIDIndex = 0,
                                      .strand = dataTypes::GenomicStrand::FORWARD,
                                      .start = 39,
                                      .end = 46,
                                      .complementarityScore = 1.0,
                                      .hybridizationEnergy = -1.7,
                                      .crosslinkingSiteCount = 1}},
     {"SRR18331301.5", "SRR18331301.4"},
     {1.0, 1.0},
     {-1.7, -1.7},
     {1, 1}},
    std::nullopt, std::nullopt);

INSTANTIATE_TEST_SUITE_P(Default, InteractionClusterGeneratorTests,
                         testing::Values(InteractionClusterGeneratorTestParam{
                             .expectedFeatureCounts = {{"gene1", 3}, {"gene3", 3}},
                             .finalizedClusters = {cluster1},
                             .partiallyAnnotatedClusters = {cluster2, cluster3}}));
