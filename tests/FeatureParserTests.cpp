// NOLINTBEGIN

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

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

    auto directGroups = FeatureGrouper::groupByDirectParentID(std::vector<GenomicFeature>{features});
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
        annotationFile
            << "##gff-version 3\n"
            << "chr1\t.\texon\t101\t150\t.\t+\t.\t"
            << "ID=exon:ENST000001.1:1;gene_id=ENSG000001.1;transcript_id=ENST000001.1\n"
            << "chr1\t.\texon\t201\t250\t.\t+\t.\t"
            << "ID=exon:ENST000001.1:2;gene_id=ENSG000001.1;transcript_id=ENST000001.1\n";
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    const auto results =
        FeatureParser{{"exon"}, std::nullopt}.parseFlatAndGrouped(annotationPath, referenceIndex);

    ASSERT_TRUE(results.flatByChromosomeIndex.contains(0));
    ASSERT_EQ(results.flatByChromosomeIndex.at(0).size(), 2);
    EXPECT_EQ(results.flatByChromosomeIndex.at(0)[0].getAnnotationID(), "ENST000001.1");
    EXPECT_EQ(results.flatByChromosomeIndex.at(0)[1].getAnnotationID(), "ENST000001.1");

    ASSERT_TRUE(results.groupedByParentID.contains("ENST000001.1"));
    EXPECT_EQ(results.groupedByParentID.at("ENST000001.1").getNodes().size(), 2);

    std::filesystem::remove(annotationPath);
}

TEST(FeatureParserTests, ExonIdentifierFallbackRestoresTranscriptGrouping) {
    const auto annotationPath =
        std::filesystem::temp_directory_path() / "rnanue_feature_parser_exon_id_fallback.gff";

    {
        std::ofstream annotationFile{annotationPath};
        annotationFile
            << "##gff-version 3\n"
            << "chr1\t.\texon\t101\t150\t.\t+\t.\tID=exon:ENST000002.1:1;gene_id=ENSG000002.1\n"
            << "chr1\t.\texon\t201\t250\t.\t+\t.\tID=exon:ENST000002.1:2;gene_id=ENSG000002.1\n";
    }

    const ReferenceIndexMapping referenceIndex{{{"chr1", 0}}, MissingReferencePolicy::Reject};
    const auto results =
        FeatureParser{{"exon"}, std::nullopt}.parseFlatAndGrouped(annotationPath, referenceIndex);

    ASSERT_TRUE(results.flatByChromosomeIndex.contains(0));
    ASSERT_EQ(results.flatByChromosomeIndex.at(0).size(), 2);
    EXPECT_EQ(results.flatByChromosomeIndex.at(0)[0].getAnnotationID(), "ENST000002.1");
    EXPECT_EQ(results.flatByChromosomeIndex.at(0)[1].getAnnotationID(), "ENST000002.1");

    ASSERT_TRUE(results.groupedByParentID.contains("ENST000002.1"));
    EXPECT_EQ(results.groupedByParentID.at("ENST000002.1").getNodes().size(), 2);

    std::filesystem::remove(annotationPath);
}

// NOLINTEND
