#include <gtest/gtest.h>

#include "gtest/gtest.h"

// Standard
#include <string>
#include <unordered_map>
#include <vector>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "EvaluatedInteractionCluster.hpp"
#include "GenomicStrand.hpp"
#include "Orientation.hpp"
#include "StatisticEvaluator.hpp"

using namespace pipelines::analyze;

// Note: These are the expected results
// Two clusters, one intra and one intermolecular
// Inter feature A:B, Intra feature C:C
// {{"A", 15UL}, {"B", 15UL}, {"C", 25UL}}

struct StatisticEvaluatorTestParam {
    std::vector<AnnotatedInteractionCluster> clusters;
    std::unordered_map<std::string, size_t> transcriptFrequencies;
    size_t totalTranscripts;
    double padjThreshold;
    std::vector<EvaluatedInteractionCluster> expectedResults;
};

class StatisticEvaluatorTests : public ::testing::TestWithParam<StatisticEvaluatorTestParam> {};

TEST_P(StatisticEvaluatorTests, Default) {
    const auto& param = GetParam();

    auto clusters = param.clusters;

    auto results = StatisticEvaluator::evaluate(clusters, param.transcriptFrequencies,
                                                param.totalTranscripts, param.padjThreshold);

    ASSERT_EQ(results.size(), param.expectedResults.size());

    for (const auto& result : results) {
        std::cout << "Result: " << result.getFirstSegmentRecordIDs().front()
                  << ", pvalue: " << result.getPValue() << ", padj: " << result.getPadj() << "\n";
    }
}

const std::unordered_map<std::string, size_t> transcriptFrequencies1{
    {"A", 15UL}, {"B", 15UL}, {"C", 25UL}};

const InteractionCluster intraMolecularCluster(
    InteractionSegment(0, Strand::FORWARD, 0, 10), InteractionSegment(0, Strand::FORWARD, 20, 30),
    {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"},
    {"11", "12", "13", "14", "15", "16", "17", "18", "19", "20"},
    {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5},
    {-200, -200, -200, -200, -200, -200, -200, -200, -200, -200});

const InteractionCluster interMolecularCluster(
    InteractionSegment(0, Strand::FORWARD, 30, 40), InteractionSegment(0, Strand::FORWARD, 45, 55),
    {"21", "22", "23", "24", "25", "26", "27", "28", "29", "30"},
    {"31", "32", "33", "34", "35", "36", "37", "38", "39", "40"},
    {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5},
    {-200, -200, -200, -200, -200, -200, -200, -200, -200, -200});

const AnnotatedInteractionCluster intraMolecularClusterAnno(intraMolecularCluster, "C", "C");
const AnnotatedInteractionCluster interMolecularClusterAnno(interMolecularCluster, "A", "B");

const std::vector<AnnotatedInteractionCluster> clusters{intraMolecularClusterAnno,
                                                        interMolecularClusterAnno};

const std::vector<EvaluatedInteractionCluster> expectedResults{
    {intraMolecularClusterAnno, 0.1818, 0.1}, {interMolecularClusterAnno, 0.1818, 0.1}};

INSTANTIATE_TEST_SUITE_P(Default, StatisticEvaluatorTests,
                         testing::Values(StatisticEvaluatorTestParam{
                             clusters, transcriptFrequencies1, 55UL, 1, expectedResults}));
