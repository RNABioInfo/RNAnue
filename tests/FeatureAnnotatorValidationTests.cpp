#include <gtest/gtest.h>

#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>
#include <unordered_set>

#include "AnnotationHierarchyError.hpp"
#include "FeatureAnnotator.hpp"
#include "FeatureParser.hpp"

namespace {

class FeatureAnnotatorValidationTests : public testing::TestWithParam<std::string> {
   protected:
    std::filesystem::path annotationPath;
    const annotation::ReferenceIDToIndexMap references{{"chr1", 0}};

    void SetUp() override {
        // Masked annotations have a .gtf suffix but are written as GFF3.
        annotationPath = std::filesystem::temp_directory_path() /
                         ("rnanue_annotator_" +
                          boost::uuids::to_string(boost::uuids::random_generator{}()) + ".gtf");
        std::ofstream out{annotationPath};
        ASSERT_TRUE(out.is_open());
        out << "##gff-version 3\n"
            << "chr1\t.\t" << GetParam() << "\t1\t100\t.\t+\t.\tID=root\n"
            << "chr1\t.\ttranscript\t1\t100\t.\t+\t.\tID=tx;Parent=root\n"
            << "chr1\t.\texon\t1\t80\t.\t+\t.\tID=ex;Parent=tx\n";
    }

    void TearDown() override {
        std::error_code error;
        std::filesystem::remove(annotationPath, error);
    }
};

TEST_P(FeatureAnnotatorValidationTests, ExcludedParentPropagatesInsteadOfTerminating) {
    const std::unordered_set<std::string> types{"transcript", "exon"};
    try {
        const annotation::FeatureAnnotator annotator{annotationPath, references, types};
        FAIL() << "Expected excluded-parent validation to fail";
    } catch (const annotation::AnnotationHierarchyError& error) {
        const std::string message{error.what()};
        EXPECT_NE(message.find("add '" + GetParam() + "' to --featuretypes"), std::string::npos);
    }
}

TEST_P(FeatureAnnotatorValidationTests, ExplicitIdConstructorAlsoPropagatesValidation) {
    const std::unordered_set<std::string> types{"transcript", "exon"};
    EXPECT_THROW((annotation::FeatureAnnotator{annotationPath, references, types, "ID"}),
                 annotation::AnnotationHierarchyError);
}

TEST_P(FeatureAnnotatorValidationTests, IncludingParentLoadsTheCompleteHierarchy) {
    const std::unordered_set<std::string> types{GetParam(), "transcript", "exon"};
    const annotation::FeatureAnnotator annotator{annotationPath, references, types};
    EXPECT_EQ(annotator.featureCount(), 3);
    const annotation::FeatureAnnotator explicitId{annotationPath, references, types, "ID"};
    EXPECT_EQ(explicitId.featureCount(), 3);
}

TEST_P(FeatureAnnotatorValidationTests, PreflightRejectsExcludedParentsWithoutABamDictionary) {
    const annotation::FeatureParser parser{{"transcript", "exon"}, std::nullopt};
    EXPECT_THROW(parser.validateHierarchy(annotationPath), annotation::AnnotationHierarchyError);
}

TEST_P(FeatureAnnotatorValidationTests, PreflightAcceptsCompleteHierarchyAcrossReferences) {
    {
        std::ofstream out{annotationPath, std::ios::app};
        out << "chr2\t.\ttranscript\t1\t100\t.\t+\t.\tID=tx2\n"
            << "chr2\t.\texon\t1\t80\t.\t+\t.\tID=ex2;Parent=tx2\n";
    }
    const annotation::FeatureParser parser{{GetParam(), "transcript", "exon"}, std::nullopt};
    EXPECT_NO_THROW(parser.validateHierarchy(annotationPath));
}

INSTANTIATE_TEST_SUITE_P(ParentTypes, FeatureAnnotatorValidationTests,
                        testing::Values("pseudogene", "primary_transcript", "scaRNA"));

}  // namespace
