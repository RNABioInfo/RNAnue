// NOLINTBEGIN

#include <gtest/gtest.h>

#include <cstddef>

// Standard
#include <iostream>
#include <optional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

// Internal
#include "FeatureAnnotator.hpp"
#include "GenomicFeature.hpp"
#include "GenomicFeatureTreeMerger.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "Region.hpp"

using namespace annotation;

struct FeatureAnnotatorTestParameters {
    FeatureAnnotatorTestParameters(dataTypes::GenomicRegion region,
                                   dataTypes::GenomicOrientation orientation,
                                   std::vector<std::string> expectedFeatureIds)
        : region(region),
          expectedFeatureIds(std::move(expectedFeatureIds)),
          orientation(orientation) {}

    FeatureAnnotatorTestParameters(dataTypes::GenomicRegion region,
                                   std::vector<std::string> expectedFeatureIds)
        : region(region), expectedFeatureIds(std::move(expectedFeatureIds)) {}

    dataTypes::GenomicRegion region;
    std::vector<std::string> expectedFeatureIds;
    dataTypes::GenomicOrientation orientation = dataTypes::GenomicOrientation::SAME;
};

// Tests for FeatureAnnotator::overlappingFeatures and FeatureAnnotator::overlappingFeatureIterator
class FeatureAnnotatorTest : public testing::TestWithParam<FeatureAnnotatorTestParameters> {
   protected:
    FeatureAnnotatorTest() : annotator(featureMap) {}

    const dataTypes::FeatureMap featureMap = {
        {
            1,
            {GenomicFeature{"transcript",
                            GenomicRegion{1, Region{.startPosition = 1, .endPosition = 10},
                                          dataTypes::GenomicStrand::FORWARD},
                            "feature1", "group1", std::nullopt},
             GenomicFeature{"transcript",
                            GenomicRegion{1, Region{.startPosition = 20, .endPosition = 30},
                                          dataTypes::GenomicStrand::FORWARD},
                            "feature2", "group1", std::nullopt},
             GenomicFeature{"transcript",
                            GenomicRegion{1, Region{.startPosition = 40, .endPosition = 50},
                                          dataTypes::GenomicStrand::FORWARD},
                            "feature3", "group2", std::nullopt},
             GenomicFeature{"transcript",
                            GenomicRegion{1, Region{.startPosition = 5, .endPosition = 25},
                                          dataTypes::GenomicStrand::REVERSE},
                            "feature4", "group3", std::nullopt}},
        },
        {2,
         {GenomicFeature{"transcript",
                         GenomicRegion{2, Region{.startPosition = 1, .endPosition = 10},
                                       dataTypes::GenomicStrand::FORWARD},
                         "feature3", "group4", std::nullopt},
          GenomicFeature{"transcript",
                         GenomicRegion{2, Region{.startPosition = 20, .endPosition = 30},
                                       dataTypes::GenomicStrand::FORWARD},
                         "feature4", "group4", std::nullopt}}},
    };

    FeatureAnnotator annotator;
};

void PrintTo(const FeatureAnnotatorTestParameters& param, std::ostream* outputStream) {
    *outputStream << param.region << "; " << param.orientation << ";\nExpected Feature IDs: ";

    for (const auto& featureID : param.expectedFeatureIds) {
        *outputStream << featureID << ", ";
    }

    *outputStream << "\n";
}

TEST_P(FeatureAnnotatorTest, OverlappingFeatures) {
    const auto& param = GetParam();
    const auto features = annotator.getOverlappingFeatures(param.region, param.orientation);

    ASSERT_EQ(features.size(), param.expectedFeatureIds.size());

    for (size_t i = 0; i < features.size(); ++i) {
        EXPECT_EQ(features[i].getID(), param.expectedFeatureIds[i]);
    }
}

TEST_P(FeatureAnnotatorTest, OverlappingFeatureIterator) {
    const auto& param = GetParam();
    const auto results = annotator.overlappingFeatureIt(param.region, param.orientation);

    size_t index = 0;
    for (const auto& feature : results) {
        EXPECT_EQ(feature.getID(), param.expectedFeatureIds[index]);
        ++index;
    }
    EXPECT_EQ(index, param.expectedFeatureIds.size());
}

INSTANTIATE_TEST_SUITE_P(
    Default, FeatureAnnotatorTest,
    testing::Values(
        FeatureAnnotatorTestParameters{
            GenomicRegion{1, Region{9, 15}, dataTypes::GenomicStrand::FORWARD}, {"feature1"}},
        FeatureAnnotatorTestParameters{
            GenomicRegion{1, Region{5, 15}, dataTypes::GenomicStrand::REVERSE}, {"feature4"}},
        FeatureAnnotatorTestParameters{GenomicRegion{1, {5, 15}, dataTypes::GenomicStrand::FORWARD},
                                       GenomicOrientation::BOTH,
                                       {"feature1", "feature4"}},
        FeatureAnnotatorTestParameters{
            GenomicRegion{1, Region{31, 39}, dataTypes::GenomicStrand::FORWARD}, {}},
        FeatureAnnotatorTestParameters{
            GenomicRegion{2, Region{5, 15}, dataTypes::GenomicStrand::FORWARD}, {"feature3"}},
        FeatureAnnotatorTestParameters{
            GenomicRegion{3, Region{5, 25}, dataTypes::GenomicStrand::FORWARD}, {}},
        FeatureAnnotatorTestParameters{
            GenomicRegion{1, Region{5, 15}, dataTypes::GenomicStrand::FORWARD},
            dataTypes::GenomicOrientation::OPPOSITE,
            {"feature4"}}));

// Tests for FeatureAnnotator::getBestOverlappingFeature
class BestFeatureAnnotatorTest : public testing::TestWithParam<FeatureAnnotatorTestParameters> {
   protected:
    BestFeatureAnnotatorTest() : annotator(featureMap) {}

    const dataTypes::FeatureMap featureMap = {
        {1,
         {GenomicFeature{"transcript",
                         GenomicRegion{1, Region{.startPosition = 1, .endPosition = 10},
                                       dataTypes::GenomicStrand::FORWARD},
                         "feature1", "group1", std::nullopt},
          GenomicFeature{"transcript",
                         GenomicRegion{1, Region{.startPosition = 20, .endPosition = 30},
                                       dataTypes::GenomicStrand::FORWARD},
                         "feature2", "group1", std::nullopt},
          GenomicFeature{"transcript",
                         GenomicRegion{1, Region{.startPosition = 40, .endPosition = 50},
                                       dataTypes::GenomicStrand::FORWARD},
                         "feature3", "group2", std::nullopt},
          GenomicFeature{"transcript",
                         GenomicRegion{1, Region{.startPosition = 5, .endPosition = 25},
                                       dataTypes::GenomicStrand::REVERSE},
                         "feature4", "group3", std::nullopt}}},
        {2,
         {GenomicFeature{"transcript",
                         GenomicRegion{2, Region{.startPosition = 1, .endPosition = 10},
                                       dataTypes::GenomicStrand::FORWARD},
                         "feature3", "group4", std::nullopt},
          GenomicFeature{"transcript",
                         GenomicRegion{2, Region{.startPosition = 20, .endPosition = 30},
                                       dataTypes::GenomicStrand::FORWARD},
                         "feature4", "group4", std::nullopt}}},
    };

    FeatureAnnotator annotator;
};

TEST_P(BestFeatureAnnotatorTest, GetBestOverlappingFeature) {
    const auto& param = GetParam();
    const auto feature =
        annotator.getBestOverlappingFeature(param.region, dataTypes::GenomicOrientation::BOTH);

    if (param.expectedFeatureIds.empty()) {
        EXPECT_FALSE(feature.has_value());
    } else {
        EXPECT_TRUE(feature.has_value());
        EXPECT_EQ(feature->getID(), param.expectedFeatureIds[0]);
    }
}

INSTANTIATE_TEST_SUITE_P(
    Default, BestFeatureAnnotatorTest,
    testing::Values(
        FeatureAnnotatorTestParameters{GenomicRegion{1, {1, 7}, dataTypes::GenomicStrand::FORWARD},
                                       {"feature1"}},
        FeatureAnnotatorTestParameters{GenomicRegion{1, {5, 15}, dataTypes::GenomicStrand::REVERSE},
                                       {"feature4"}},
        FeatureAnnotatorTestParameters{GenomicRegion{1, {5, 15}, dataTypes::GenomicStrand::NONE},
                                       {"feature4"}},
        FeatureAnnotatorTestParameters{GenomicRegion{1, {31, 39}, dataTypes::GenomicStrand::NONE},
                                       {}},
        FeatureAnnotatorTestParameters{GenomicRegion{2, {5, 15}, dataTypes::GenomicStrand::NONE},
                                       {"feature3"}},
        FeatureAnnotatorTestParameters{GenomicRegion{3, {5, 25}, dataTypes::GenomicStrand::NONE},
                                       {}}));

// Tests for FeatureAnnotator::insert and FeatureAnnotator::mergeInsert
class InsertFeatureAnnotatorTest : public testing::Test {
   protected:
    InsertFeatureAnnotatorTest() : annotator(featureMap) {}

    const dataTypes::FeatureMap featureMap = {
        {1,
         {
             GenomicFeature{"transcript",
                            GenomicRegion{1, Region{.startPosition = 1, .endPosition = 10},
                                          dataTypes::GenomicStrand::FORWARD},
                            "feature1", "group1", std::nullopt},
             GenomicFeature{"transcript",
                            GenomicRegion{1, Region{.startPosition = 20, .endPosition = 30},
                                          dataTypes::GenomicStrand::FORWARD},
                            "feature2", "group1", std::nullopt},
         }},
    };

    FeatureAnnotator annotator;
};

TEST_F(InsertFeatureAnnotatorTest, Insert) {
    const GenomicRegion region{
        1, {.startPosition = 5, .endPosition = 19}, dataTypes::GenomicStrand::FORWARD};
    annotator.insertIndex(region);

    const auto features =
        annotator.getOverlappingFeatures(region, dataTypes::GenomicOrientation::SAME);
    ASSERT_EQ(features.size(), 2UL);
    EXPECT_EQ(features[0].getID(), "feature1");
    EXPECT_NE(features[1].getID(), "feature2");
}

class MergeFeatureAnnotatorTest : public testing::Test {
   protected:
    MergeFeatureAnnotatorTest() : annotator(featureMap) {}

    const dataTypes::FeatureMap featureMap = {
        {1,
         {GenomicFeature{
              "transcript",
              {1, {.startPosition = 1, .endPosition = 10}, dataTypes::GenomicStrand::FORWARD},
              "feature1",
              "group1",
              std::nullopt},
          GenomicFeature{
              "transcript",
              {1, {.startPosition = 20, .endPosition = 30}, dataTypes::GenomicStrand::FORWARD},
              "feature2",
              "group1",
              std::nullopt},
          GenomicFeature{
              "transcript",
              {1, {.startPosition = 8, .endPosition = 50}, dataTypes::GenomicStrand::FORWARD},
              "feature3",
              "group2",
              std::nullopt},
          GenomicFeature{
              "transcript",
              {1, {.startPosition = 5, .endPosition = 25}, dataTypes::GenomicStrand::REVERSE},
              "feature4",
              "group3",
              std::nullopt},
          GenomicFeature{
              "transcript",
              {1, {.startPosition = 50, .endPosition = 56}, dataTypes::GenomicStrand::FORWARD},
              "feature5",
              "group2",
              std::nullopt}}},
        {2,
         {GenomicFeature{
              "transcript",
              {2, {.startPosition = 1, .endPosition = 10}, dataTypes::GenomicStrand::FORWARD},
              "feature6",
              "group4",
              std::nullopt},
          GenomicFeature{
              "transcript",
              {2, {.startPosition = 8, .endPosition = 30}, dataTypes::GenomicStrand::FORWARD},
              "feature7",
              "group4",
              std::nullopt}}},
    };

    FeatureAnnotator annotator;
};

TEST_F(MergeFeatureAnnotatorTest, MergeOverlapOne) {
    FeatureMergingParameters parameters{GenomicStrandSpecificity::SPECIFIC, -1};
    annotator.mergeAllOverlappingFeatures(parameters);

    ASSERT_EQ(annotator.featureCount(), 4UL);

    const auto regionOne = GenomicRegion{
        1, {.startPosition = 1, .endPosition = 50}, dataTypes::GenomicStrand::FORWARD};
    const auto feature1 = annotator.getOverlappingFeatures(regionOne, GenomicOrientation::SAME);

    ASSERT_EQ(feature1.size(), 1UL);
    EXPECT_EQ(feature1[0].getID(), "feature1");
    EXPECT_EQ(feature1[0].getGenomicRegion().getStart(), 1);
    EXPECT_EQ(feature1[0].getGenomicRegion().getEnd(), 50);

    const auto regionTwo = GenomicRegion{
        1, {.startPosition = 5, .endPosition = 25}, dataTypes::GenomicStrand::REVERSE};
    const auto feature2 = annotator.getOverlappingFeatures(regionTwo, GenomicOrientation::SAME);
    ASSERT_EQ(feature2.size(), 1UL);
    EXPECT_EQ(feature2[0].getID(), "feature4");
    EXPECT_EQ(feature2[0].getGenomicRegion().getStart(), 5);
    EXPECT_EQ(feature2[0].getGenomicRegion().getEnd(), 25);

    const auto regionThree = GenomicRegion{
        1, {.startPosition = 51, .endPosition = 56}, dataTypes::GenomicStrand::FORWARD};
    const auto feature3 = annotator.getOverlappingFeatures(regionThree, GenomicOrientation::SAME);
    ASSERT_EQ(feature3.size(), 1UL);
    EXPECT_EQ(feature3[0].getID(), "feature5");
    EXPECT_EQ(feature3[0].getGenomicRegion().getStart(), 50);
    EXPECT_EQ(feature3[0].getGenomicRegion().getEnd(), 56);

    const auto regionFour = GenomicRegion{
        2, {.startPosition = 1, .endPosition = 30}, dataTypes::GenomicStrand::FORWARD};
    const auto feature4 = annotator.getOverlappingFeatures(regionFour, GenomicOrientation::SAME);
    ASSERT_EQ(feature4.size(), 1UL);
    EXPECT_EQ(feature4[0].getID(), "feature6");
    EXPECT_EQ(feature4[0].getGenomicRegion().getStart(), 1);
    EXPECT_EQ(feature4[0].getGenomicRegion().getEnd(), 30);
}

// NOLINTEND
