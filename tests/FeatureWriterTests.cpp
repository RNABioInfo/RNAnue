#include <gtest/gtest.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "FeatureParser.hpp"
#include "FeatureWriter.hpp"
#include "GenomicFeatureGroup.hpp"
#include "GenomicStrand.hpp"
#include "MaskedFeatureCluster.hpp"

using namespace annotation;
using namespace dataTypes;

namespace {
auto feature(std::string id, std::optional<std::string> parent = std::nullopt)
    -> GenomicFeature {
    return {"transcript", GenomicRegion{0, {0, 10}, GenomicStrand::FORWARD},
            std::move(id), std::move(parent), std::nullopt};
}

class FeatureWriterTests : public ::testing::Test {
protected:
    std::filesystem::path directory;
    std::filesystem::path output;

    void SetUp() override {
        directory = std::filesystem::temp_directory_path() /
                    ("rnanue-writer-" + std::to_string(
                        std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directory(directory);
        output = directory / "masked.gff";
    }

    void TearDown() override { std::filesystem::remove_all(directory); }

    auto readOutput() const -> std::string {
        std::ifstream input{output};
        return {std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    }
};
}

TEST_F(FeatureWriterTests, OrphanForestFailsBeforeTruncatingOutput) {
    // Same structure as the two-locus RNAcentral groups in the server crash.
    auto group = GenomicFeatureGroup::buildFromFlat(
        {feature("rnacentral-locus-1", "URS0001BD83D9"),
         feature("rnacentral-locus-2", "URS0001BD83D9")}, "URS0001BD83D9");
    ASSERT_EQ(group.group.getNodes().size(), 2);
    ASSERT_EQ(group.group.getNodes()[1].parentIndex, GenomicFeatureGroup::invalidIndex);
    std::vector<MaskedFeatureCluster> clusters{{std::move(group.group), {}, {}}};
    { std::ofstream existing{output}; existing << "existing output\n"; }
    EXPECT_THROW(FeatureWriter::write(clusters, {"chr1"}, output.string(), FileType::GFF),
                 std::invalid_argument);
    EXPECT_EQ(readOutput(), "existing output\n");
}

TEST_F(FeatureWriterTests, ValidHierarchyRoundTripsCanonicalParentsAndMaskedCopies) {
    auto base = GenomicFeatureGroup::buildFromFlat(
        {feature("root"), feature("child", "root"), feature("grandchild", "child")});
    auto copy = GenomicFeatureGroup::buildFromFlat({feature("copy")});
    std::vector<MaskedFeatureCluster> clusters{
        {std::move(base.group), {}, {std::move(copy.group)}}};
    FeatureWriter::write(clusters, {"chr1"}, output.string(), FileType::GFF);
    EXPECT_NE(readOutput().find("MaskedFeatures=copy"), std::string::npos);
    const ReferenceIndexMapping references{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    const auto groups = FeatureParser::defaultAll().parseGroupedByHierarchy(output, references);
    ASSERT_EQ(groups.size(), 1);
    const auto& group = groups.at("root");
    ASSERT_EQ(group.getNodes().size(), 3);
    const auto& grandchild = group.getNode(*group.tryFindById("grandchild"));
    EXPECT_EQ(group.getNode(grandchild.parentIndex).feature.getID(), "child");
}

TEST_F(FeatureWriterTests, InvalidMaskedCopyFailsBeforeCreatingOutput) {
    auto base = GenomicFeatureGroup::buildFromFlat({feature("root")});
    std::vector<MaskedFeatureCluster> clusters{{std::move(base.group), {}, {GenomicFeatureGroup{}}}};
    EXPECT_THROW(FeatureWriter::write(clusters, {"chr1"}, output.string(), FileType::GFF),
                 std::invalid_argument);
    EXPECT_FALSE(std::filesystem::exists(output));
}

TEST_F(FeatureWriterTests, SharedRnaCentralMetadataKeepsLociSeparate) {
    const auto inputPath = directory / "rnacentral.gff";
    {
        std::ofstream input{inputPath};
        input << "##gff-version 3\n"
              << "chr1\tRNAcentral\tsnRNA\t1\t10\t.\t+\t.\t"
                 "ID=rnacentral-locus-1;Name=URS0001BD83D9;aliases=URS0001BD83D9;"
                 "transcript_id=URS0001BD83D9;gene_id=URS0001BD83D9\n"
              << "chr1\tRNAcentral\tsnRNA\t21\t30\t.\t-\t.\t"
                 "ID=rnacentral-locus-2;Name=URS0001BD83D9;aliases=URS0001BD83D9;"
                 "transcript_id=URS0001BD83D9;gene_id=URS0001BD83D9\n";
    }
    const ReferenceIndexMapping references{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    auto groups = FeatureParser::defaultAllTranscript().parseGroupedByHierarchy(inputPath, references);
    ASSERT_EQ(groups.size(), 2);
    std::vector<MaskedFeatureCluster> clusters;
    for (auto& [id, group] : groups) {
        EXPECT_FALSE(group.getRoot().feature.getParentID());
        clusters.push_back({std::move(group), {}, {}});
    }
    EXPECT_NO_THROW(FeatureWriter::write(clusters, {"chr1"}, output.string(), FileType::GFF));
    const auto reparsed = FeatureParser::defaultAllTranscript().parseGroupedByHierarchy(output, references);
    EXPECT_EQ(reparsed.size(), 2);
}

TEST_F(FeatureWriterTests, CustomLongIdAttributeSurvivesBothParserEntryPoints) {
    const std::string key = "custom_feature_identifier_longer_than_small_string_storage";
    const auto inputPath = directory / "custom.gff";
    {
        std::ofstream input{inputPath};
        input << "##gff-version 3\nchr1\t.\ttranscript\t1\t10\t.\t+\t.\t"
              << key << "=root\nchr1\t.\texon\t1\t5\t.\t+\t.\t"
              << key << "=child;Parent=root\n";
    }
    const ReferenceIndexMapping references{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    const FeatureParser parser{{"transcript", "exon"}, key};
    const auto groups = parser.parseGroupedByHierarchy(inputPath, references);
    ASSERT_EQ(groups.size(), 1);
    EXPECT_EQ(groups.at("root").getNodes().size(), 2);
    const auto flat = parser.parseFlatAndGrouped(inputPath, references);
    ASSERT_EQ(flat.flatByChromosomeIndex.at(0).size(), 2);
    EXPECT_EQ(flat.flatByChromosomeIndex.at(0)[1].getID(), "child");
}

TEST(GenomicFeatureGroupSafetyTests, InvalidAccessIsCatchable) {
    const GenomicFeatureGroup empty;
    EXPECT_THROW(empty.getRoot(), std::out_of_range);
    EXPECT_THROW(empty.getNode(GenomicFeatureGroup::invalidIndex), std::out_of_range);
    EXPECT_THROW(empty.children(GenomicFeatureGroup::invalidIndex), std::out_of_range);
    EXPECT_THROW(empty.allChildrenOfRoot(), std::out_of_range);
    EXPECT_THROW(empty.requireRootedTree("test"), std::invalid_argument);
}

TEST(GenomicFeatureGroupSafetyTests, DisconnectedCycleIsRejected) {
    const auto built = GenomicFeatureGroup::buildFromFlat(
        {feature("root"), feature("a", "b"), feature("b", "a")}, "root");
    EXPECT_THROW(built.group.requireRootedTree("test"), std::invalid_argument);
}
