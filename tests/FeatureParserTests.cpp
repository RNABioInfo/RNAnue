// NOLINTBEGIN

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include "AnnotationHierarchyError.hpp"
#include "FeatureGrouper.hpp"
#include "FeatureParser.hpp"
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "ReferenceIndexMapping.hpp"
#include "Region.hpp"

using namespace annotation;
using namespace dataTypes;

namespace {

auto testRegion(int start, int end) -> GenomicRegion {
    return GenomicRegion{0, Region{.startPosition = start, .endPosition = end},
                         GenomicStrand::FORWARD};
}

}  // namespace

TEST(FeatureGrouperTests, DirectParentGroupingDoesNotCollapseToHierarchyRoot) {
    std::vector<GenomicFeature> features{
        GenomicFeature{"gene", testRegion(0, 1), "gene1", std::nullopt, std::nullopt},
        GenomicFeature{"transcript", testRegion(0, 10), "tx1", "gene1", std::nullopt},
        GenomicFeature{"exon", testRegion(0, 4), "exon:tx1:1", "tx1", std::nullopt},
        GenomicFeature{"exon", testRegion(6, 10), "exon:tx1:2", "tx1", std::nullopt},
    };

    auto directGroups =
        FeatureGrouper::groupByDirectParentID(std::vector<GenomicFeature>{features});
    auto hierarchyGroups = FeatureGrouper::groupByHierarchy(std::move(features));

    EXPECT_TRUE(directGroups.groups.contains("tx1"));
    EXPECT_EQ(directGroups.groups.at("tx1").getNodes().size(), 2);
    EXPECT_EQ(directGroups.groups.size(), 2);

    EXPECT_TRUE(hierarchyGroups.groups.contains("gene1"));
    EXPECT_EQ(hierarchyGroups.groups.size(), 1);
}

TEST(FeatureParserTests, ExonRecordsWithoutParentFallBackToTranscriptId) {
    const auto annotationPath =
        std::filesystem::temp_directory_path() / "rnanue_feature_parser_parent_fallback.gff";

    {
        std::ofstream annotationFile{annotationPath};
        annotationFile << "##gff-version 3\n"
                       << "chr1\t.\ttranscript\t101\t250\t.\t+\t.\tID=ENST000001.1\n"
                       << "chr1\t.\texon\t101\t150\t.\t+\t.\t"
                       << "ID=exon:ENST000001.1:1;gene_id=ENSG000001.1;transcript_id="
                          "ENST000001.1\n"
                       << "chr1\t.\texon\t201\t250\t.\t+\t.\t"
                       << "ID=exon:ENST000001.1:2;gene_id=ENSG000001.1;transcript_id="
                          "ENST000001.1\n";
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    const auto results = FeatureParser{{"transcript", "exon"}, std::nullopt}.parseFlatAndGrouped(
        annotationPath, referenceIndex);

    ASSERT_TRUE(results.flatByChromosomeIndex.contains(0));
    ASSERT_EQ(results.flatByChromosomeIndex.at(0).size(), 3);
    EXPECT_EQ(results.flatByChromosomeIndex.at(0)[1].getAnnotationID(), "ENST000001.1");
    EXPECT_EQ(results.flatByChromosomeIndex.at(0)[2].getAnnotationID(), "ENST000001.1");

    ASSERT_TRUE(results.groupedByParentID.contains("ENST000001.1"));
    EXPECT_EQ(results.groupedByParentID.at("ENST000001.1").getNodes().size(), 3);

    std::filesystem::remove(annotationPath);
}

TEST(FeatureParserTests, ExonIdentifierFallbackRestoresTranscriptGrouping) {
    const auto annotationPath =
        std::filesystem::temp_directory_path() / "rnanue_feature_parser_exon_id_fallback.gff";

    {
        std::ofstream annotationFile{annotationPath};
        annotationFile << "##gff-version 3\n"
                       << "chr1\t.\ttranscript\t101\t250\t.\t+\t.\tID=ENST000002.1\n"
                       << "chr1\t.\texon\t101\t150\t.\t+\t.\tID=exon:ENST000002.1:1;gene_id="
                          "ENSG000002.1\n"
                       << "chr1\t.\texon\t201\t250\t.\t+\t.\tID=exon:ENST000002.1:2;gene_id="
                          "ENSG000002.1\n";
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    const auto results = FeatureParser{{"transcript", "exon"}, std::nullopt}.parseFlatAndGrouped(
        annotationPath, referenceIndex);

    ASSERT_TRUE(results.flatByChromosomeIndex.contains(0));
    ASSERT_EQ(results.flatByChromosomeIndex.at(0).size(), 3);
    EXPECT_EQ(results.flatByChromosomeIndex.at(0)[1].getAnnotationID(), "ENST000002.1");
    EXPECT_EQ(results.flatByChromosomeIndex.at(0)[2].getAnnotationID(), "ENST000002.1");

    ASSERT_TRUE(results.groupedByParentID.contains("ENST000002.1"));
    EXPECT_EQ(results.groupedByParentID.at("ENST000002.1").getNodes().size(), 3);

    std::filesystem::remove(annotationPath);
}

TEST(FeatureParserTests, NativeTranscriptMetadataDoesNotBecomeParent) {
    const auto annotationPath =
        std::filesystem::temp_directory_path() / "rnanue_native_transcript_root.gff";
    {
        std::ofstream annotationFile{annotationPath};
        annotationFile << "##gff-version 3\n"
                       << "chr1\tGENCODE\ttranscript\t101\t250\t.\t+\t.\t"
                       << "ID=gencode-ENST000001.1;transcript_id=ENST000001.1;gene_id="
                          "ENSG000001.1\n"
                       << "chr1\tGENCODE\texon\t101\t150\t.\t+\t.\t"
                       << "ID=gencode-ENST000001.1-exon-1;Parent=gencode-ENST000001.1\n";
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    const auto groups = FeatureParser{{"transcript", "exon"}, std::nullopt}.parseGroupedByHierarchy(
        annotationPath, referenceIndex);

    ASSERT_TRUE(groups.contains("gencode-ENST000001.1"));
    ASSERT_EQ(groups.size(), 1);
    const auto& group = groups.at("gencode-ENST000001.1");
    EXPECT_EQ(group.getRoot().feature.getID(), "gencode-ENST000001.1");
    EXPECT_FALSE(group.getRoot().feature.getParentID().has_value());
    EXPECT_EQ(group.getNodes().size(), 2);

    std::filesystem::remove(annotationPath);
}

TEST(FeatureParserTests, ExplicitParentResolvesThroughUniqueDeclaredAlias) {
    const auto annotationPath =
        std::filesystem::temp_directory_path() / "rnanue_parent_alias_repair.gff";
    {
        std::ofstream annotationFile{annotationPath};
        annotationFile << "##gff-version 3\n"
                       << "chr1\tGENCODE\ttranscript\t101\t250\t.\t+\t.\t"
                       << "ID=gencode-ENST000002.1;aliases=ENST000002.1\n"
                       << "chr1\tGENCODE\texon\t101\t150\t.\t+\t.\t"
                       << "ID=exon-1;Parent=ENST000002.1\n";
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    const auto groups = FeatureParser{{"transcript", "exon"}, std::nullopt}.parseGroupedByHierarchy(
        annotationPath, referenceIndex);

    ASSERT_TRUE(groups.contains("gencode-ENST000002.1"));
    ASSERT_EQ(groups.at("gencode-ENST000002.1").getNodes().size(), 2);

    std::filesystem::remove(annotationPath);
}

TEST(FeatureParserTests, MissingExonParentIsInferredFromStructuredIdentifier) {
    const auto annotationPath =
        std::filesystem::temp_directory_path() / "rnanue_structured_exon_parent.gff";
    {
        std::ofstream annotationFile{annotationPath};
        annotationFile << "##gff-version 3\n"
                       << "chr1\t.\ttranscript\t101\t250\t.\t+\t.\tID=ENST000003.1\n"
                       << "chr1\t.\texon\t101\t150\t.\t+\t.\tID=exon:ENST000003.1:1\n";
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    const auto groups = FeatureParser{{"transcript", "exon"}, std::nullopt}.parseGroupedByHierarchy(
        annotationPath, referenceIndex);

    ASSERT_TRUE(groups.contains("ENST000003.1"));
    EXPECT_EQ(groups.at("ENST000003.1").getNodes().size(), 2);

    std::filesystem::remove(annotationPath);
}

TEST(FeatureParserTests, ExcludedParentFailsWithFeatureTypeGuidance) {
    const auto annotationPath =
        std::filesystem::temp_directory_path() / "rnanue_excluded_parent.gff";
    {
        std::ofstream annotationFile{annotationPath};
        annotationFile << "##gff-version 3\n"
                       << "chr1\t.\tgene\t101\t250\t.\t+\t.\tID=gene-1\n"
                       << "chr1\t.\texon\t101\t150\t.\t+\t.\tID=exon-1;Parent=gene-1\n";
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    try {
        static_cast<void>(FeatureParser{{"exon"}, std::nullopt}.parseGroupedByHierarchy(
            annotationPath, referenceIndex));
        FAIL() << "Expected annotation hierarchy validation to fail";
    } catch (const AnnotationHierarchyError& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("add 'gene' to --featuretypes"), std::string::npos);
        EXPECT_NE(message.find("Line 3"), std::string::npos);
    }

    std::filesystem::remove(annotationPath);
}

TEST(FeatureParserTests, AmbiguousAliasFailsDescriptively) {
    const auto annotationPath =
        std::filesystem::temp_directory_path() / "rnanue_ambiguous_parent_alias.gff";
    {
        std::ofstream annotationFile{annotationPath};
        annotationFile << "##gff-version 3\n"
                       << "chr1\t.\ttranscript\t101\t250\t.\t+\t.\tID=tx-1;aliases=shared\n"
                       << "chr1\t.\ttranscript\t101\t250\t.\t+\t.\tID=tx-2;aliases=shared\n"
                       << "chr1\t.\texon\t101\t150\t.\t+\t.\tID=exon-1;Parent=shared\n";
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    try {
        static_cast<void>(
            FeatureParser{{"transcript", "exon"}, std::nullopt}.parseGroupedByHierarchy(
                annotationPath, referenceIndex));
        FAIL() << "Expected ambiguous alias validation to fail";
    } catch (const AnnotationHierarchyError& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("ambiguous between 'tx-1'"), std::string::npos);
        EXPECT_NE(message.find("'tx-2'"), std::string::npos);
    }

    std::filesystem::remove(annotationPath);
}

TEST(FeatureParserTests, InvalidHierarchyAggregatesMultipleFailures) {
    const auto annotationPath =
        std::filesystem::temp_directory_path() / "rnanue_invalid_hierarchy.gff";
    {
        std::ofstream annotationFile{annotationPath};
        annotationFile << "##gff-version 3\n"
                       << "chr1\t.\ttranscript\t101\t250\t.\t+\t.\tID=cycle-a;Parent=cycle-b\n"
                       << "chr1\t.\ttranscript\t101\t250\t.\t+\t.\tID=cycle-b;Parent=cycle-a\n"
                       << "chr1\t.\texon\t101\t150\t.\t+\t.\tID=multi;Parent=cycle-a,cycle-b\n"
                       << "chr1\t.\texon\t251\t260\t.\t+\t.\tID=outside;Parent=cycle-a\n";
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    try {
        static_cast<void>(
            FeatureParser{{"transcript", "exon"}, std::nullopt}.parseGroupedByHierarchy(
                annotationPath, referenceIndex));
        FAIL() << "Expected invalid hierarchy validation to fail";
    } catch (const AnnotationHierarchyError& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("multiple Parents"), std::string::npos);
        EXPECT_NE(message.find("outside Parent"), std::string::npos);
        EXPECT_NE(message.find("cycle detected"), std::string::npos);
        EXPECT_NE(message.find("No masking or alignment was performed"), std::string::npos);
    }

    std::filesystem::remove(annotationPath);
}

TEST(FeatureParserTests, EquivalentDuplicateRecordsAreDeduplicated) {
    const auto annotationPath =
        std::filesystem::temp_directory_path() / "rnanue_equivalent_duplicate.gff";
    {
        std::ofstream annotationFile{annotationPath};
        annotationFile << "##gff-version 3\n"
                       << "chr1\t.\ttranscript\t101\t250\t.\t+\t.\tID=tx-duplicate\n"
                       << "chr1\t.\ttranscript\t101\t250\t.\t+\t.\tID=tx-duplicate\n";
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    const auto groups = FeatureParser{{"transcript"}, std::nullopt}.parseGroupedByHierarchy(
        annotationPath, referenceIndex);
    ASSERT_TRUE(groups.contains("tx-duplicate"));
    EXPECT_EQ(groups.at("tx-duplicate").getNodes().size(), 1);

    std::filesystem::remove(annotationPath);
}

TEST(FeatureParserTests, GtfTranscriptGroupingRemainsUnchanged) {
    const auto annotationPath = std::filesystem::temp_directory_path() / "rnanue_gtf_grouping.gtf";
    {
        std::ofstream annotationFile{annotationPath};
        annotationFile << "##gtf-version 2.2\n"
                       << "chr1\t.\texon\t101\t150\t.\t+\t.\t"
                       << "gene_id \"gene-1\"; transcript_id \"tx-1\"; gene \"GENE1\";\n";
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    const auto groups = FeatureParser{{"exon"}, std::nullopt}.parseGroupedByParentID(
        annotationPath, referenceIndex);
    ASSERT_TRUE(groups.contains("tx-1"));
    EXPECT_EQ(groups.at("tx-1").getNodes().size(), 1);

    std::filesystem::remove(annotationPath);
}

TEST(FeatureParserTests, LargeExactHierarchyScalesWithoutRepairIndex) {
    const auto annotationPath =
        std::filesystem::temp_directory_path() / "rnanue_large_exact_hierarchy.gff";
    constexpr std::size_t transcriptCount = 5'000;
    {
        std::ofstream annotationFile{annotationPath};
        annotationFile << "##gff-version 3\n";
        for (std::size_t index = 0; index < transcriptCount; ++index) {
            const auto start = index * 3 + 1;
            const auto end = start + 1;
            annotationFile << "chr1\t.\ttranscript\t" << start << '\t' << end
                           << "\t.\t+\t.\tID=tx-" << index << '\n'
                           << "chr1\t.\texon\t" << start << '\t' << end
                           << "\t.\t+\t.\tID=exon-" << index << ";Parent=tx-" << index
                           << '\n';
        }
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    const auto groups = FeatureParser{{"transcript", "exon"}, std::nullopt}
                            .parseGroupedByHierarchy(annotationPath, referenceIndex);
    EXPECT_EQ(groups.size(), transcriptCount);

    std::filesystem::remove(annotationPath);
}

TEST(FeatureParserTests, ErrorReportBoundsStoredDiagnostics) {
    const auto annotationPath =
        std::filesystem::temp_directory_path() / "rnanue_bounded_hierarchy_errors.gff";
    constexpr std::size_t errorCount = 60;
    {
        std::ofstream annotationFile{annotationPath};
        annotationFile << "##gff-version 3\n";
        for (std::size_t index = 0; index < errorCount; ++index) {
            annotationFile << "chr1\t.\texon\t1\t1\t.\t+\t.\tID=exon-" << index
                           << ";Parent=missing-" << index << '\n';
        }
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    try {
        static_cast<void>(FeatureParser{{"exon"}, std::nullopt}.parseGroupedByHierarchy(
            annotationPath, referenceIndex));
        FAIL() << "Expected bounded hierarchy diagnostics";
    } catch (const AnnotationHierarchyError& error) {
        const std::string message = error.what();
        EXPECT_NE(message.find("with 60 errors"), std::string::npos);
        EXPECT_NE(message.find("10 additional errors omitted"), std::string::npos);
    }

    std::filesystem::remove(annotationPath);
}

// NOLINTEND
