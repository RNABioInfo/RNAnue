// NOLINTBEGIN(readability-magic-numbers)

#include <gtest/gtest.h>

// Standard
#include <deque>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "EvaluatedInteractionCluster.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "InteractionCluster.hpp"
#include "InteractionsWriter.hpp"
#include "RecordFragment.hpp"
#include "Region.hpp"

using namespace dataTypes;
using namespace pipelines::analyze;

namespace {

auto readFile(const std::filesystem::path& path) -> std::string {
    std::ifstream in{path};
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

auto makeFragment(std::string recordID, int referenceID, int start, int end)
    -> RecordFragment {
    return {.recordID = std::move(recordID),
            .genomicRegion =
                GenomicRegion{referenceID, Region{.startPosition = start, .endPosition = end},
                              GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -1.0,
            .interCrosslinkingSiteCount = 1,
            .transcriptContribution = 1.0F,
            .coverageIntervals = {{.start = start, .end = end}}};
}

auto makeEvaluatedCluster() -> EvaluatedInteractionCluster {
    auto cluster = InteractionCluster::fromRecordFragments(makeFragment("read1", 0, 0, 10),
                                                           makeFragment("read1", 1, 100, 110));
    return EvaluatedInteractionCluster{
        AnnotatedInteractionCluster{std::move(cluster), "featureA", "featureB"}, 0.01, 0.02};
}

}  // namespace

TEST(InteractionsWriterTests, WritesCoverageMetricsAndAggregateBedGraph) {
    const auto testDir =
        std::filesystem::temp_directory_path() / "RNAnueInteractionsWriterTests";
    std::filesystem::remove_all(testDir);
    std::filesystem::create_directories(testDir);

    const InteractionsWriter::OutputPaths outputPaths{
        .interactionsOutputPath = testDir / "sample_interactions.tsv",
        .interactionReadIDsOutputPath = testDir / "sample_interactions_readIDs.tsv",
        .interactionsBEDOutputPath = testDir / "sample_interaction_regions.bed",
        .interactionsBEDArcOutputPath = testDir / "sample_interaction_regions.arc",
        .interactionArmCoverageBedGraphOutputPath =
            testDir / "sample_interaction_arm_coverage.bedgraph"};

    const std::deque<std::string> referenceIDs{"chr1", "chr2"};
    const std::vector<EvaluatedInteractionCluster> clusters{makeEvaluatedCluster()};

    InteractionsWriter::writeInteractions("sample", outputPaths, referenceIDs, clusters);

    const std::string interactions = readFile(outputPaths.interactionsOutputPath);
    EXPECT_NE(interactions.find("support_per_effective_bp"), std::string::npos);
    EXPECT_NE(interactions.find("coverage_profile"), std::string::npos);
    EXPECT_NE(interactions.find("compact_dense"), std::string::npos);

    const std::string bedGraph = readFile(outputPaths.interactionArmCoverageBedGraphOutputPath);
    EXPECT_NE(bedGraph.find("track type=bedGraph"), std::string::npos);
    EXPECT_NE(bedGraph.find("chr1\t0\t10\t1.000000"), std::string::npos);
    EXPECT_NE(bedGraph.find("chr2\t100\t110\t1.000000"), std::string::npos);

    std::filesystem::remove_all(testDir);
}

// NOLINTEND(readability-magic-numbers)
