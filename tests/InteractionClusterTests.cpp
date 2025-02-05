// NOLINTBEGIN(readability-magic-numbers)

#include <gtest/gtest.h>

// Internal
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "InteractionCluster.hpp"
#include "RecordFragment.hpp"
#include "Region.hpp"

using namespace pipelines::analyze;
using namespace dataTypes;

class InteractionClusterStatisticsTest : public testing::Test {
   protected:
    InteractionClusterStatisticsTest()
        : cluster(InteractionCluster::fromRecordFragments(
              RecordFragment{
                  .recordID = "record1",
                  .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                                 GenomicStrand::FORWARD},
                  .complementarityScore = 0.5,
                  .hybridizationEnergy = -18,
                  .interCrosslinkingSiteCount = 1},
              RecordFragment{
                  .recordID = "record1",
                  .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                                 GenomicStrand::FORWARD},
                  .complementarityScore = 0.5,
                  .hybridizationEnergy = -18,
                  .interCrosslinkingSiteCount = 1})) {};

    InteractionCluster cluster;
};

TEST_F(InteractionClusterStatisticsTest, ComplementarityStatistics) {
    EXPECT_DOUBLE_EQ(cluster.complementarityStatistics(), 0.25);
}

TEST_F(InteractionClusterStatisticsTest, HybridizationEnergyStatistics) {
    EXPECT_DOUBLE_EQ(cluster.hybridizationEnergyStatistics(), 18.0);
}

// InteractionCluster overlap tests
TEST(InteractionClusterTest, OverlapsExact) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    EXPECT_TRUE(cluster1.overlapsWithTolerance(cluster2, GenomicStrandSpecificity::SPECIFIC, -9));
}

TEST(InteractionClusterTest, OverlapsGrace) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 10, .endPosition = 20},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    InteractionCluster cluster3 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 11, .endPosition = 20},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    EXPECT_TRUE(cluster1.overlapsWithTolerance(cluster2, GenomicStrandSpecificity::SPECIFIC, 1));
    EXPECT_FALSE(cluster1.overlapsWithTolerance(cluster3, GenomicStrandSpecificity::SPECIFIC, 1));
}

TEST(InteractionClusterTest, OverlapsNegativeGrace) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 8, .endPosition = 20},
                                           GenomicStrand::REVERSE},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::REVERSE},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    EXPECT_TRUE(cluster1.overlapsWithTolerance(cluster2, GenomicStrandSpecificity::UNSPECIFIC, -1));
}

TEST(InteractionClusterTest, NoOverlapsNegativeGrace) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::REVERSE},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::REVERSE},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 9, .endPosition = 20},
                                           GenomicStrand::REVERSE},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::REVERSE},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    EXPECT_FALSE(cluster1.overlapsWithTolerance(cluster2, GenomicStrandSpecificity::SPECIFIC, -2));
}

TEST(InteractionClusterTest, NoOverlapsFirstNotOverlapping) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});
    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 11, .endPosition = 20},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    EXPECT_FALSE(cluster1.overlapsWithTolerance(cluster2, GenomicStrandSpecificity::SPECIFIC, -1));
}

TEST(InteractionClusterTest, NoOverlapsSecondNotOverlapping) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 10, .endPosition = 20},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    EXPECT_FALSE(cluster1.overlapsWithTolerance(cluster2, GenomicStrandSpecificity::SPECIFIC, 0));
}

TEST(InteractionClusterTest, NoOverlapsDifferentReferenceID) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{2, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    EXPECT_FALSE(cluster1.overlapsWithTolerance(cluster2, GenomicStrandSpecificity::UNSPECIFIC, 0));
}

// InteractionCluster merge tests
TEST(InteractionClusterTest, BasicMerge) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 5, .endPosition = 15},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 20, .endPosition = 25},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});

    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.8,
            .hybridizationEnergy = -15,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 20, .endPosition = 30},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.8,
            .hybridizationEnergy = -15,
            .interCrosslinkingSiteCount = 1});

    EXPECT_TRUE(cluster1.merge(cluster2, GenomicStrandSpecificity::SPECIFIC));

    EXPECT_EQ(cluster1.fragmentCount(), 2UL);
    EXPECT_EQ(cluster1.getComplementarityScores().size(), 2UL);
    EXPECT_EQ(cluster1.getHybridizationEnergies().size(), 2UL);

    EXPECT_EQ(cluster1.getFirstSegment().getStart(), 0);
    EXPECT_EQ(cluster1.getFirstSegment().getEnd(), 15);

    EXPECT_EQ(cluster1.getSecondSegment().getStart(), 20);
    EXPECT_EQ(cluster1.getSecondSegment().getEnd(), 30);

    EXPECT_EQ(cluster1.getFirstSegment().getReferenceIDIndex(), 0);
    EXPECT_EQ(cluster1.getSecondSegment().getReferenceIDIndex(), 0);

    EXPECT_EQ(cluster1.getMinHybridizationEnergy(), -18);

    EXPECT_EQ(cluster1.getMaxComplementarityScore(), 0.8);
}

TEST(InteractionClusterTest, MergeWithDifferentReferenceIDIndex) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 5, .endPosition = 15},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 20, .endPosition = 25},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});
    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 20, .endPosition = 30},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.8,
            .hybridizationEnergy = -15,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{1, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.8,
            .hybridizationEnergy = -15,
            .interCrosslinkingSiteCount = 1});

    EXPECT_TRUE(cluster1.merge(cluster2, GenomicStrandSpecificity::SPECIFIC));

    EXPECT_EQ(cluster1.fragmentCount(), 2UL);
    EXPECT_EQ(cluster1.getComplementarityScores().size(), 2UL);
    EXPECT_EQ(cluster1.getHybridizationEnergies().size(), 2UL);

    EXPECT_EQ(cluster1.getFirstSegment().getStart(), 20);
    EXPECT_EQ(cluster1.getFirstSegment().getEnd(), 30);

    EXPECT_EQ(cluster1.getSecondSegment().getStart(), 0);
    EXPECT_EQ(cluster1.getSecondSegment().getEnd(), 15);

    EXPECT_EQ(cluster1.getFirstSegment().getReferenceIDIndex(), 0);
    EXPECT_EQ(cluster1.getSecondSegment().getReferenceIDIndex(), 1);

    EXPECT_EQ(cluster1.getMinHybridizationEnergy(), -18);

    EXPECT_EQ(cluster1.getMaxComplementarityScore(), 0.8);

    EXPECT_EQ(cluster1.getFirstSegment().getStrand(), GenomicStrand::FORWARD);
    EXPECT_EQ(cluster1.getSecondSegment().getStrand(), GenomicStrand::FORWARD);
}

TEST(InteractionClusterTest, MergeWithDifferentStrands) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 5, .endPosition = 15},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record1",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 20, .endPosition = 25},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -18,
            .interCrosslinkingSiteCount = 1});
    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 0, .endPosition = 10},
                                           GenomicStrand::REVERSE},
            .complementarityScore = 0.8,
            .hybridizationEnergy = -15,
            .interCrosslinkingSiteCount = 1},
        RecordFragment{
            .recordID = "record2",
            .genomicRegion = GenomicRegion{0, Region{.startPosition = 20, .endPosition = 30},
                                           GenomicStrand::FORWARD},
            .complementarityScore = 0.8,
            .hybridizationEnergy = -15,
            .interCrosslinkingSiteCount = 1});

    EXPECT_FALSE(cluster1.merge(cluster2, GenomicStrandSpecificity::SPECIFIC));
}

// NOLINTEND(readability-magic-numbers)
