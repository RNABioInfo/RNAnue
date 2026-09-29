#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <numeric>
#include <random>
#include <tuple>

#include "InteractionClusterComponentBuilder.hpp"
#include "InteractionsWriter.hpp"
#include "pipelines/postprocess/InteractionParser.hpp"
#include "pipelines/postprocess/PostprocessData.hpp"

using namespace dataTypes;
using namespace pipelines::analyze;

namespace {
auto observation(int i, int start = 0, GenomicStrand strand = GenomicStrand::FORWARD)
    -> InteractionCluster {
    return {SortedGenomicRegionPair{GenomicRegion{0, {start, start + 40}, strand},
                                    GenomicRegion{1, {1000 + start, 1040 + start}, strand}},
            "read" + std::to_string(i),
            0.1 * (i % 8),
            -double(i + 1),
            i % 5,
            0.25 * (1 + i % 4)};
}
using Evidence = std::tuple<std::string, double, double, int, double>;
auto evidence(const InteractionCluster& cluster) -> std::vector<Evidence> {
    std::vector<Evidence> out;
    for (size_t i = 0; i < cluster.fragmentCount(); ++i)
        out.emplace_back(cluster.getRecordIDs()[i], cluster.getComplementarityScores()[i],
                         cluster.getHybridizationEnergies()[i],
                         cluster.getCrosslinkingSiteCounts()[i], cluster.getContributions()[i]);
    std::sort(out.begin(), out.end());
    return out;
}
auto exhaustive(std::vector<InteractionCluster> input, const ClusteringParameters& params)
    -> std::vector<InteractionCluster> {
    // Independent, deliberately slow graph traversal; no interval index or compression.
    std::vector<bool> visited(input.size());
    std::vector<InteractionCluster> result;
    for (size_t seed = 0; seed < input.size(); ++seed) {
        if (visited[seed]) continue;
        std::vector<size_t> members{seed};
        visited[seed] = true;
        for (size_t cursor = 0; cursor < members.size(); ++cursor) {
            for (size_t candidate = 0; candidate < input.size(); ++candidate) {
                if (!visited[candidate] &&
                    (params.clusterMergingStrandSpecificity ==
                         GenomicStrandSpecificity::UNSPECIFIC ||
                     input[members[cursor]].getFirstSegment().getStrand() ==
                         input[candidate].getFirstSegment().getStrand()) &&
                    InteractionClusterComponentBuilder::clustersOverlap(input[members[cursor]],
                                                                        input[candidate], params)) {
                    visited[candidate] = true;
                    members.push_back(candidate);
                }
            }
        }
        auto cluster = input[seed];
        for (size_t i = 1; i < members.size(); ++i)
            cluster.absorbValidatedComponentMember(input[members[i]]);
        result.push_back(std::move(cluster));
    }
    return result;
}
auto byMembership(const std::vector<InteractionCluster>& clusters) {
    std::map<std::vector<Evidence>, const InteractionCluster*> out;
    for (const auto& cluster : clusters) out.emplace(evidence(cluster), &cluster);
    return out;
}
}  // namespace

TEST(ClusteringRegression, MatchesExhaustiveComponentsAndPreservesEveryMetric) {
    std::vector<InteractionCluster> input;
    for (int i = 0; i < 30; ++i)
        input.push_back(
            observation(i, 30 * (i % 7), i % 4 ? GenomicStrand::FORWARD : GenomicStrand::REVERSE));
    for (auto strand : {GenomicStrandSpecificity::SPECIFIC, GenomicStrandSpecificity::UNSPECIFIC}) {
        for (ClusteringMergeParameterVariant parameter :
             std::vector<ClusteringMergeParameterVariant>{
                 ClusterOverlapToleranceMergeParameter{1},
                 ClusterOverlapToleranceMergeParameter{-50},
                 ShortestClusterOverlapFractionMergeParameter{0.3F},
                 ShortestClusterOverlapFractionMergeParameter{1.0F}}) {
            ClusteringParameters params{parameter, strand, 1, 0, GenomicOrientation::BOTH};
            auto expected = exhaustive(input, params);
            std::shuffle(input.begin(), input.end(), std::mt19937{20260924});
            std::vector<InteractionCluster> actual;
            for (auto& group : InteractionClusterComponentBuilder::groupClusters(
                     std::vector<InteractionCluster>(input), strand)) {
                auto result =
                    InteractionClusterComponentBuilder::mergeGroup(std::move(group), params);
                for (auto& cluster : result) actual.push_back(std::move(cluster));
            }
            const auto expectedMap = byMembership(expected), actualMap = byMembership(actual);
            ASSERT_EQ(actualMap.size(), expectedMap.size());
            for (const auto& [members, cluster] : expectedMap) {
                auto found = actualMap.find(members);
                ASSERT_NE(found, actualMap.end());
                const auto& other = *found->second;
                EXPECT_EQ(cluster->getFirstSegment(), other.getFirstSegment());
                EXPECT_EQ(cluster->getSecondSegment(), other.getSecondSegment());
                EXPECT_EQ(cluster->getFirstArmCoverageRuns(), other.getFirstArmCoverageRuns());
                EXPECT_EQ(cluster->getSecondArmCoverageRuns(), other.getSecondArmCoverageRuns());
                EXPECT_DOUBLE_EQ(cluster->getTranscriptContribution(),
                                 other.getTranscriptContribution());
                EXPECT_DOUBLE_EQ(cluster->coverageShapeMetrics().effectiveCoverageSpanBp,
                                 other.coverageShapeMetrics().effectiveCoverageSpanBp);
                EXPECT_EQ(cluster->coverageShapeMetrics().coverageProfile,
                          other.coverageShapeMetrics().coverageProfile);
            }
        }
    }
}

TEST(ClusteringRegression, DenseDuplicatesNeedLinearOverlapComparisons) {
    std::vector<InteractionCluster> input;
    for (int i = 0; i < 4000; ++i) input.push_back(observation(i));
    ClusteringParameters params{ClusterOverlapToleranceMergeParameter{1},
                                GenomicStrandSpecificity::SPECIFIC, 1, 0, GenomicOrientation::BOTH};
    InteractionClusterComponentBuilder::ComparisonMetrics metrics;
    auto clusters =
        InteractionClusterComponentBuilder::mergeGroup(std::move(input), params, &metrics);
    ASSERT_EQ(clusters.size(), 1);
    EXPECT_EQ(clusters.front().fragmentCount(), 4000);
    EXPECT_EQ(metrics.indexedCount, 1);
    EXPECT_EQ(metrics.overlapComparisons, 3999);
}

TEST(PostprocessingRegression, ScoresAndSmallNumbersRoundTripWithoutSchemaChange) {
    const auto dir =
        std::filesystem::temp_directory_path() / ("RNAnue-roundtrip-" + helper::getUUID());
    std::filesystem::create_directories(dir);
    auto cluster = InteractionCluster{
        SortedGenomicRegionPair{GenomicRegion{0, {0, 40}, GenomicStrand::FORWARD},
                                GenomicRegion{1, {100, 140}, GenomicStrand::FORWARD}},
        "one",
        0.5,
        -7,
        1,
        1.0 / 3};
    const auto gcs = cluster.complementarityStatistics(),
               ghs = cluster.hybridizationEnergyStatistics();
    std::vector<EvaluatedInteractionCluster> clusters;
    clusters.emplace_back(AnnotatedInteractionCluster{std::move(cluster), "A", "B"}, 1e-100,
                          2e-100);
    InteractionsWriter::OutputPaths paths{dir / "sample_interactions.tsv", dir / "ids.tsv",
                                          dir / "regions.bed", dir / "arcs.bed",
                                          dir / "coverage.bedgraph"};
    InteractionsWriter::writeInteractions("sample", paths, {"chr1", "chr2"}, clusters);
    pipelines::postprocess::InteractionParser parser;
    parser.parse(paths.interactionsOutputPath, "sample");
    auto result = parser.getResultAndReset();
    ASSERT_EQ(result.interactions.size(), 1);
    const auto& metrics = result.interactions.front().getInteractionMetrics().front();
    EXPECT_DOUBLE_EQ(metrics.globalComplementarityScore, gcs);
    EXPECT_DOUBLE_EQ(metrics.globalHybridizationScore, ghs);
    EXPECT_DOUBLE_EQ(metrics.contributionScore, 1.0 / 3);
    EXPECT_NEAR(metrics.pValue, 1e-100, 1e-112);
    EXPECT_NEAR(metrics.padjValue, 2e-100, 1e-112);
    std::filesystem::remove_all(dir);
}

TEST(PostprocessingRegression, DuplicateSampleNamesFailAcrossConditions) {
    const auto dir =
        std::filesystem::temp_directory_path() / ("RNAnue-samples-" + helper::getUUID());
    for (const auto& condition : {"treatment", "control"}) {
        std::filesystem::create_directories(dir / condition / "same");
        std::ofstream(dir / condition / "same/same_interactions.tsv") << "unused\n";
    }
    EXPECT_THROW(
        (pipelines::postprocess::PostprocessData{dir / "out", dir / "treatment", dir / "control"}),
        std::runtime_error);
    std::filesystem::remove_all(dir);
}
