#include <gtest/gtest.h>

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "AnnotatedInteractionCluster.hpp"
#include "FeatureAnnotator.hpp"
#include "GenomicStrand.hpp"
#include "InteractionCluster.hpp"
#include "InteractionClusterGenerator.hpp"
#include "Orientation.hpp"
#include "RecordFragment.hpp"
#include "SplitRecordsParser.hpp"

using namespace pipelines::analyze;

// Note: This is the expected clustering result:
// Cluster 1: 1,3,6
// Cluster 2: 2,7
// Cluster 3: 4,5

using namespace pipelines::analyze;

auto testInteractionClustersSamPath() -> std::string {
    return (std::filesystem::path{__FILE__}.parent_path() /
            "test_data/interactionClusterRecords.sam")
        .string();
}

struct InteractionClusterGeneratorTestParam {
    std::vector<InteractionCluster> expectedInteractionClusters;
};

std::shared_ptr<annotation::FeatureAnnotator> featureAnnotator =
    std::make_shared<annotation::FeatureAnnotator>();

class InteractionClusterGeneratorTests
    : public testing::TestWithParam<InteractionClusterGeneratorTestParam> {
   protected:
    InteractionClusterGeneratorTests()
        : interactionClusterGenerator("Testing Sample", featureAnnotator, {"chr1", "chr2"},
                                      annotation::Orientation::BOTH, 0, 1) {}

    InteractionClusterGenerator interactionClusterGenerator;
};

// TODO: Fix this test and use == operator for InteractionCluster2
TEST_P(InteractionClusterGeneratorTests, SplitRecordsAreSortedCorrectly) {
    const auto& param = GetParam();
    const auto& expectedClusters = param.expectedInteractionClusters;

    std::vector<InteractionCluster> interactionClusters =
        SplitRecordsParser::parse(testInteractionClustersSamPath());

    EXPECT_EQ(interactionClusters.size(), 7U);

    auto result = interactionClusterGenerator.mergeClusters(std::move(interactionClusters));

    EXPECT_EQ(result.annotatedClusters.size(), expectedClusters.size());

    for (const auto& mergedCluster : result.annotatedClusters) {
        std::cout << "Merged cluster: " << mergedCluster << "\n";

        bool found = false;

        for (const auto& expectedCluster : expectedClusters) {
            if (mergedCluster.getFirstSegment().getStart() ==
                    expectedCluster.getFirstSegment().getStart() &&
                mergedCluster.getFirstSegment().getEnd() ==
                    expectedCluster.getFirstSegment().getEnd() &&
                mergedCluster.getSecondSegment().getStart() ==
                    expectedCluster.getSecondSegment().getStart() &&
                mergedCluster.getSecondSegment().getEnd() ==
                    expectedCluster.getSecondSegment().getEnd() &&
                mergedCluster.getFirstSegment().getReferenceIDIndex() ==
                    expectedCluster.getFirstSegment().getReferenceIDIndex() &&
                mergedCluster.getSecondSegment().getReferenceIDIndex() ==
                    expectedCluster.getSecondSegment().getReferenceIDIndex() &&
                mergedCluster.getFirstSegment().getStrand() ==
                    expectedCluster.getFirstSegment().getStrand() &&
                mergedCluster.getSecondSegment().getStrand() ==
                    expectedCluster.getSecondSegment().getStrand()) {
                found = true;
                break;
            }
        }

        EXPECT_TRUE(found);
    }
}

const InteractionCluster cluster1 = {RecordFragment{.recordID = "SRR18331301.3",
                                                    .referenceIDIndex = 0,
                                                    .strand = dataTypes::Strand::FORWARD,
                                                    .start = 19,
                                                    .end = 27,
                                                    .complementarityScore = 1,
                                                    .hybridizationEnergy = -1.7},
                                     RecordFragment{.recordID = "SRR18331301.3",
                                                    .referenceIDIndex = 1,
                                                    .strand = dataTypes::Strand::FORWARD,
                                                    .start = 49,
                                                    .end = 64,
                                                    .complementarityScore = 1,
                                                    .hybridizationEnergy = -1.7},
                                     {"SRR18331301.20", "SRR18331301.21", "SRR18331301.22"},
                                     {"SRR18331301.20", "SRR18331301.21", "SRR18331301.22"},
                                     {1, 1, 1},
                                     {-1.7, -1.7, -1.7}};

const InteractionCluster cluster2 = {RecordFragment{.recordID = "SRR18331301.2",
                                                    .referenceIDIndex = 1,
                                                    .strand = dataTypes::Strand::FORWARD,
                                                    .start = 4,
                                                    .end = 10,
                                                    .complementarityScore = 1,
                                                    .hybridizationEnergy = -1.7},
                                     RecordFragment{.recordID = "SRR18331301.2",
                                                    .referenceIDIndex = 1,
                                                    .strand = dataTypes::Strand::FORWARD,
                                                    .start = 51,
                                                    .end = 57,
                                                    .complementarityScore = 1,
                                                    .hybridizationEnergy = -1.7},
                                     "SRR18331301.2",
                                     "SRR18331301.2",
                                     1,
                                     -1.7};

const InteractionCluster cluster3 = {RecordFragment{.recordID = "SRR18331301.4",
                                                    .referenceIDIndex = 0,
                                                    .strand = dataTypes::Strand::FORWARD,
                                                    .start = 4,
                                                    .end = 14,
                                                    .complementarityScore = 1,
                                                    .hybridizationEnergy = -1.7},
                                     RecordFragment{.recordID = "SRR18331301.4",
                                                    .referenceIDIndex = 0,
                                                    .strand = dataTypes::Strand::FORWARD,
                                                    .start = 39,
                                                    .end = 46,
                                                    .complementarityScore = 1,
                                                    .hybridizationEnergy = -1.7},
                                     "SRR18331301.4",
                                     "SRR18331301.4",
                                     1,
                                     -1.7};

INSTANTIATE_TEST_SUITE_P(Default, InteractionClusterGeneratorTests,
                         testing::Values(InteractionClusterGeneratorTestParam{
                             {cluster1, cluster2, cluster3}}));
