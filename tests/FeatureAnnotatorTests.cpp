#include <gtest/gtest.h>

// Standard
#include <optional>
#include <utility>

// Internal
#include "FeatureAnnotator.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "Orientation.hpp"

using namespace annotation;

struct FeatureAnnotatorTestParameters {
    FeatureAnnotatorTestParameters(dataTypes::GenomicRegion region,
                                   annotation::Orientation orientation,
                                   std::vector<std::string> expectedFeatureIds)
        : region(std::move(region)),
          expectedFeatureIds(std::move(expectedFeatureIds)),
          orientation(orientation) {}

    FeatureAnnotatorTestParameters(dataTypes::GenomicRegion region,
                                   std::vector<std::string> expectedFeatureIds)
        : region(std::move(region)), expectedFeatureIds(std::move(expectedFeatureIds)) {}

    dataTypes::GenomicRegion region;
    std::vector<std::string> expectedFeatureIds;
    annotation::Orientation orientation = annotation::Orientation::SAME;
};

// Tests for FeatureAnnotator::overlappingFeatures and FeatureAnnotator::overlappingFeatureIterator
class FeatureAnnotatorTest : public testing::TestWithParam<FeatureAnnotatorTestParameters> {
   protected:
    FeatureAnnotatorTest() : annotator(featureMap) {}

    const dataTypes::FeatureMap featureMap = {
        {
            "chromosome1",
            {{.referenceID = "chromosome1",
              .type = "transcript",
              .startPosition = 1,
              .endPosition = 10,
              .strand = dataTypes::GenomicStrand::FORWARD,
              .id = "feature1",
              .groupID = "group1",
              .geneName = std::nullopt},
             {.referenceID = "chromosome1",
              .type = "transcript",
              .startPosition = 20,
              .endPosition = 30,
              .strand = dataTypes::GenomicStrand::FORWARD,
              .id = "feature2",
              .groupID = "group1",
              .geneName = std::nullopt},
             {.referenceID = "chromosome1",
              .type = "transcript",
              .startPosition = 40,
              .endPosition = 50,
              .strand = dataTypes::GenomicStrand::FORWARD,
              .id = "feature3",
              .groupID = "group2",
              .geneName = std::nullopt},
             {.referenceID = "chromosome1",
              .type = "transcript",
              .startPosition = 5,
              .endPosition = 25,
              .strand = dataTypes::GenomicStrand::REVERSE,
              .id = "feature4",
              .groupID = "group3",
              .geneName = std::nullopt}},
        },
        {"chromosome2",
         {{.referenceID = "chromosome2",
           .type = "transcript",
           .startPosition = 1,
           .endPosition = 10,
           .strand = dataTypes::GenomicStrand::FORWARD,
           .id = "feature3",
           .groupID = "group4",
           .geneName = std::nullopt},
          {.referenceID = "chromosome2",
           .type = "transcript",
           .startPosition = 20,
           .endPosition = 30,
           .strand = dataTypes::GenomicStrand::FORWARD,
           .id = "feature4",
           .groupID = "group4",
           .geneName = std::nullopt}}},
    };

    FeatureAnnotator annotator;
};

void PrintTo(const FeatureAnnotatorTestParameters& param, std::ostream* outputStream) {
    *outputStream << "GenomicRegion{referenceID: " << param.region.referenceID
                  << ", startPosition: " << param.region.startPosition
                  << ", endPosition: " << param.region.endPosition << ", strand: "
                  << (param.region.strand.has_value() ? std::to_string(param.region.strand.value())
                                                      : "nullopt")
                  << "}";
}

TEST_P(FeatureAnnotatorTest, OverlappingFeatures) {
    const auto& param = GetParam();
    const auto features = annotator.getOverlappingFeatures(param.region, param.orientation);

    ASSERT_EQ(features.size(), param.expectedFeatureIds.size());
    for (size_t i = 0; i < features.size(); ++i) {
        EXPECT_EQ(features[i].id, param.expectedFeatureIds[i]);
    }
}

TEST_P(FeatureAnnotatorTest, OverlappingFeatureIterator) {
    const auto& param = GetParam();
    const auto results = annotator.overlappingFeatureIterator(param.region, param.orientation);

    size_t i = 0;
    for (const auto& feature : results) {
        EXPECT_EQ(feature.id, param.expectedFeatureIds[i]);
        ++i;
    }
    EXPECT_EQ(i, param.expectedFeatureIds.size());
}

INSTANTIATE_TEST_SUITE_P(
    Default, FeatureAnnotatorTest,
    testing::Values(
        FeatureAnnotatorTestParameters{{"chromosome1", 9, 15, dataTypes::GenomicStrand::FORWARD},
                                       {"feature1"}},
        FeatureAnnotatorTestParameters{{"chromosome1", 5, 15, dataTypes::GenomicStrand::REVERSE},
                                       {"feature4"}},
        FeatureAnnotatorTestParameters{{"chromosome1", 5, 15, std::nullopt},
                                       {"feature1", "feature4"}},
        FeatureAnnotatorTestParameters{{"chromosome1", 31, 39, std::nullopt}, {}},
        FeatureAnnotatorTestParameters{{"chromosome2", 5, 15, std::nullopt}, {"feature3"}},
        FeatureAnnotatorTestParameters{{"chromosome3", 5, 25, std::nullopt}, {}},
        FeatureAnnotatorTestParameters{{"chromosome1", 5, 15, dataTypes::GenomicStrand::FORWARD},
                                       annotation::Orientation::OPPOSITE,
                                       {"feature4"}}));

// Tests for FeatureAnnotator::getBestOverlappingFeature
class BestFeatureAnnotatorTest : public testing::TestWithParam<FeatureAnnotatorTestParameters> {
   protected:
    BestFeatureAnnotatorTest() : annotator(featureMap) {}

    const dataTypes::FeatureMap featureMap = {
        {"chromosome1",
         {{"chromosome1", "transcript", 1, 10, dataTypes::GenomicStrand::FORWARD, "feature1",
           "group1", std::nullopt},
          {"chromosome1", "transcript", 20, 30, dataTypes::GenomicStrand::FORWARD, "feature2",
           "group1", std::nullopt},
          {"chromosome1", "transcript", 40, 50, dataTypes::GenomicStrand::FORWARD, "feature3",
           "group2", std::nullopt},
          {"chromosome1", "transcript", 5, 25, dataTypes::GenomicStrand::REVERSE, "feature4",
           "group3", std::nullopt}}},
        {"chromosome2",
         {{"chromosome2", "transcript", 1, 10, dataTypes::GenomicStrand::FORWARD, "feature3",
           "group4", std::nullopt},
          {"chromosome2", "transcript", 20, 30, dataTypes::GenomicStrand::FORWARD, "feature4",
           "group4", std::nullopt}}},
    };

    FeatureAnnotator annotator;
};

TEST_P(BestFeatureAnnotatorTest, GetBestOverlappingFeature) {
    const auto& param = GetParam();
    const auto feature =
        annotator.getBestOverlappingFeature(param.region, annotation::Orientation::BOTH);

    if (param.expectedFeatureIds.empty()) {
        EXPECT_FALSE(feature.has_value());
    } else {
        EXPECT_TRUE(feature.has_value());
        EXPECT_EQ(feature->id, param.expectedFeatureIds[0]);
    }
}

INSTANTIATE_TEST_SUITE_P(
    Default, BestFeatureAnnotatorTest,
    testing::Values(
        FeatureAnnotatorTestParameters{{"chromosome1", 1, 7, dataTypes::GenomicStrand::FORWARD},
                                       {"feature1"}},
        FeatureAnnotatorTestParameters{{"chromosome1", 5, 15, dataTypes::GenomicStrand::REVERSE},
                                       {"feature4"}},
        FeatureAnnotatorTestParameters{{"chromosome1", 5, 15, std::nullopt}, {"feature4"}},
        FeatureAnnotatorTestParameters{{"chromosome1", 31, 39, std::nullopt}, {}},
        FeatureAnnotatorTestParameters{{"chromosome2", 5, 15, std::nullopt}, {"feature3"}},
        FeatureAnnotatorTestParameters{{"chromosome3", 5, 25, std::nullopt}, {}}));

// Tests for FeatureAnnotator::insert and FeatureAnnotator::mergeInsert
class InsertFeatureAnnotatorTest : public testing::Test {
   protected:
    InsertFeatureAnnotatorTest() : annotator(featureMap) {}

    const dataTypes::FeatureMap featureMap = {
        {"chromosome1",
         {
             {"chromosome1", "transcript", 1, 10, dataTypes::GenomicStrand::FORWARD, "feature1",
              "group1", std::nullopt},
             {"chromosome1", "transcript", 20, 30, dataTypes::GenomicStrand::FORWARD, "feature2",
              "group1", std::nullopt},
         }},
    };

    FeatureAnnotator annotator;
};

TEST_F(InsertFeatureAnnotatorTest, Insert) {
    const dataTypes::GenomicRegion region{"chromosome1", 5, 20, dataTypes::GenomicStrand::FORWARD};
    const auto featureId = annotator.insertIndex(region);

    const auto features = annotator.getOverlappingFeatures(region, annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 2ul);
    EXPECT_EQ(features[0].id, "feature1");
    EXPECT_NE(features[1].id, "feature2");
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsert) {
    const dataTypes::GenomicRegion region{"chromosome1", 5, 15, dataTypes::GenomicStrand::FORWARD};
    const auto result = annotator.mergeInsertIndex(region, 0);

    ASSERT_EQ(annotator.featureCount(), 2ul);

    EXPECT_EQ(result.featureID, "feature1");
    EXPECT_TRUE(result.mergedFeatureIDs.empty());

    const auto features = annotator.getOverlappingFeatures(region, annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1ul);
    EXPECT_EQ(features[0].id, "feature1");
    EXPECT_EQ(features[0].startPosition, 1);
    EXPECT_EQ(features[0].endPosition, 15);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertTwoOverlapping) {
    const dataTypes::GenomicRegion region{"chromosome1", 5, 25, dataTypes::GenomicStrand::FORWARD};
    const auto result = annotator.mergeInsertIndex(region, 0);

    ASSERT_EQ(annotator.featureCount(), 1ul);

    EXPECT_EQ(result.featureID, "feature1");
    ASSERT_EQ(result.mergedFeatureIDs.size(), 1ul);
    EXPECT_EQ(result.mergedFeatureIDs[0], "feature2");

    const auto features = annotator.getOverlappingFeatures(region, annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1ul);
    EXPECT_EQ(features[0].id, "feature1");
    EXPECT_EQ(features[0].startPosition, 1);
    EXPECT_EQ(features[0].endPosition, 30);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertWithReverseStrand) {
    const dataTypes::GenomicRegion region{"chromosome1", 5, 15, dataTypes::GenomicStrand::REVERSE};
    const auto result = annotator.mergeInsertIndex(region, 0);

    ASSERT_EQ(annotator.featureCount(), 3ul);

    const auto features = annotator.getOverlappingFeatures(region, annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1ul);
    EXPECT_EQ(features[0].id, result.featureID);
    EXPECT_TRUE(result.mergedFeatureIDs.empty());
    EXPECT_EQ(features[0].startPosition, 5);
    EXPECT_EQ(features[0].endPosition, 15);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertWithNoStrand) {
    const dataTypes::GenomicRegion region{"chromosome1", 5, 15, std::nullopt};
    ASSERT_DEATH(annotator.mergeInsertIndex(region, 0), "Strand must be specified for insertion");
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertNotExistingReferenceID) {
    const dataTypes::GenomicRegion region{"chromosome2", 5, 25, dataTypes::GenomicStrand::FORWARD};
    const auto result = annotator.mergeInsertIndex(region, 0);

    ASSERT_EQ(annotator.featureCount(), 3ul);

    const auto features = annotator.getOverlappingFeatures(region, annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1ul);
    EXPECT_EQ(features[0].id, result.featureID);
    EXPECT_TRUE(result.mergedFeatureIDs.empty());
    EXPECT_EQ(features[0].startPosition, 5);
    EXPECT_EQ(features[0].endPosition, 25);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertGraceDistance) {
    const dataTypes::GenomicRegion region{"chromosome1", 5, 15, dataTypes::GenomicStrand::FORWARD};
    const auto result = annotator.mergeInsertIndex(region, 5);

    ASSERT_EQ(annotator.featureCount(), 1ul);

    const auto features = annotator.getOverlappingFeatures(region, annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1ul);
    EXPECT_EQ(features[0].id, result.featureID);
    EXPECT_EQ(result.mergedFeatureIDs, std::vector<std::string>{"feature2"});
    EXPECT_EQ(features[0].startPosition, 1);
    EXPECT_EQ(features[0].endPosition, 30);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertGraceDistanceNotSecondOverlapping) {
    const dataTypes::GenomicRegion region{"chromosome1", 1, 14, dataTypes::GenomicStrand::FORWARD};
    const auto result = annotator.mergeInsertIndex(region, 5);

    ASSERT_EQ(annotator.featureCount(), 2ul);

    const auto features = annotator.getOverlappingFeatures(region, annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1ul);
    EXPECT_EQ(features[0].id, result.featureID);
    EXPECT_TRUE(result.mergedFeatureIDs.empty());
    EXPECT_EQ(features[0].startPosition, 1);
    EXPECT_EQ(features[0].endPosition, 14);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertGraceDistanceOneSpace) {
    const dataTypes::GenomicRegion region{"chromosome1", 11, 15, dataTypes::GenomicStrand::FORWARD};
    const auto result = annotator.mergeInsertIndex(region, 0);

    ASSERT_EQ(annotator.featureCount(), 3ul);

    const auto features = annotator.getOverlappingFeatures(region, annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1ul);
    EXPECT_EQ(features[0].id, result.featureID);
    EXPECT_TRUE(result.mergedFeatureIDs.empty());
    EXPECT_EQ(features[0].startPosition, 11);
    EXPECT_EQ(features[0].endPosition, 15);
}

TEST_F(InsertFeatureAnnotatorTest, MergeInsertGraceDistanceBluntEnds) {
    const dataTypes::GenomicRegion region{"chromosome1", 10, 15, dataTypes::GenomicStrand::FORWARD};
    const auto result = annotator.mergeInsertIndex(region, 0);

    ASSERT_EQ(annotator.featureCount(), 2ul);

    const auto features = annotator.getOverlappingFeatures(region, annotation::Orientation::SAME);
    ASSERT_EQ(features.size(), 1ul);
    EXPECT_EQ(features[0].id, result.featureID);
    EXPECT_TRUE(result.mergedFeatureIDs.empty());
    EXPECT_EQ(features[0].startPosition, 1);
    EXPECT_EQ(features[0].endPosition, 15);
}

class MergeFeatureAnnotatorTest : public testing::Test {
   protected:
    MergeFeatureAnnotatorTest() : annotator(featureMap) {}

    const dataTypes::FeatureMap featureMap = {
        {"chromosome1",
         {{"chromosome1", "transcript", 1, 10, dataTypes::GenomicStrand::FORWARD, "feature1",
           "group1", std::nullopt},
          {"chromosome1", "transcript", 20, 30, dataTypes::GenomicStrand::FORWARD, "feature2",
           "group1", std::nullopt},
          {"chromosome1", "transcript", 8, 50, dataTypes::GenomicStrand::FORWARD, "feature3",
           "group2", std::nullopt},
          {"chromosome1", "transcript", 5, 25, dataTypes::GenomicStrand::REVERSE, "feature4",
           "group3", std::nullopt},
          {"chromosome1", "transcript", 51, 56, dataTypes::GenomicStrand::FORWARD, "feature5",
           "group2", std::nullopt}}},
        {"chromosome2",
         {{"chromosome2", "transcript", 1, 10, dataTypes::GenomicStrand::FORWARD, "feature6",
           "group4", std::nullopt},
          {"chromosome2", "transcript", 10, 30, dataTypes::GenomicStrand::FORWARD, "feature7",
           "group4", std::nullopt}}},
    };

    FeatureAnnotator annotator;
};

TEST_F(MergeFeatureAnnotatorTest, MergeOverlapOne) {
    constexpr int minOverlap = 1;
    annotator.mergeIndexAllOverlappingFeatures(minOverlap);

    for (const auto& map : annotator.getFeatureTreeMap()) {
        for (const auto& feature : map.second.intervals()) {
            std::cout << feature.data.id << " " << feature.data.startPosition << " "
                      << feature.data.endPosition << '\n';
        }
    }

    ASSERT_EQ(annotator.featureCount(), 4UL);

    const auto regionOne =
        dataTypes::GenomicRegion{"chromosome1", 1, 50, dataTypes::GenomicStrand::FORWARD};
    const auto feature1 = annotator.getOverlappingFeatures(regionOne, Orientation::SAME);
    ASSERT_EQ(feature1.size(), 1UL);
    EXPECT_EQ(feature1[0].id, "feature1");
    EXPECT_EQ(feature1[0].startPosition, 1);
    EXPECT_EQ(feature1[0].endPosition, 50);

    const auto regionTwo =
        dataTypes::GenomicRegion{"chromosome1", 5, 25, dataTypes::GenomicStrand::REVERSE};
    const auto feature2 = annotator.getOverlappingFeatures(regionTwo, Orientation::SAME);
    ASSERT_EQ(feature2.size(), 1UL);
    EXPECT_EQ(feature2[0].id, "feature4");
    EXPECT_EQ(feature2[0].startPosition, 5);
    EXPECT_EQ(feature2[0].endPosition, 25);

    const auto regionThree =
        dataTypes::GenomicRegion{"chromosome1", 51, 56, dataTypes::GenomicStrand::FORWARD};
    const auto feature3 = annotator.getOverlappingFeatures(regionThree, Orientation::SAME);
    ASSERT_EQ(feature3.size(), 1UL);
    EXPECT_EQ(feature3[0].id, "feature5");
    EXPECT_EQ(feature3[0].startPosition, 51);
    EXPECT_EQ(feature3[0].endPosition, 56);

    const auto regionFour =
        dataTypes::GenomicRegion{"chromosome2", 1, 30, dataTypes::GenomicStrand::FORWARD};
    const auto feature4 = annotator.getOverlappingFeatures(regionFour, Orientation::SAME);
    ASSERT_EQ(feature4.size(), 1UL);
    EXPECT_EQ(feature4[0].id, "feature6");
    EXPECT_EQ(feature4[0].startPosition, 1);
    EXPECT_EQ(feature4[0].endPosition, 30);
}
