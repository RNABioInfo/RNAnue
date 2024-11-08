#include <gtest/gtest.h>

// Internal
#include "GenomicStrand.hpp"
#include "InteractionCluster.hpp"

using namespace pipelines::analyze;
using namespace dataTypes;

class InteractionClusterStatisticsTest : public testing::Test {
   protected:
    InteractionClusterStatisticsTest()
        : cluster(InteractionCluster::fromRecordFragments(
              RecordFragment{.recordID = "record1",
                             .referenceIDIndex = 0,
                             .strand = GenomicStrand::FORWARD,
                             .start = 0,
                             .end = 10,
                             .complementarityScore = 0.5,
                             .hybridizationEnergy = -18},
              RecordFragment{.recordID = "record1",
                             .referenceIDIndex = 0,
                             .strand = GenomicStrand::FORWARD,
                             .start = 0,
                             .end = 10,
                             .complementarityScore = 0.5,
                             .hybridizationEnergy = -18})) {};

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
    InteractionCluster cluster1 =
        InteractionCluster::fromRecordFragments(RecordFragment{.recordID = "record1",
                                                               .referenceIDIndex = 0,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 0,
                                                               .end = 10,
                                                               .complementarityScore = 0.5,
                                                               .hybridizationEnergy = -18},
                                                RecordFragment{.recordID = "record1",
                                                               .referenceIDIndex = 1,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 0,
                                                               .end = 10,
                                                               .complementarityScore = 0.5,
                                                               .hybridizationEnergy = -18});
    InteractionCluster cluster2 =
        InteractionCluster::fromRecordFragments(RecordFragment{.recordID = "record2",
                                                               .referenceIDIndex = 0,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 0,
                                                               .end = 10,
                                                               .complementarityScore = 0.5,
                                                               .hybridizationEnergy = -18},
                                                RecordFragment{.recordID = "record2",
                                                               .referenceIDIndex = 1,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 0,
                                                               .end = 10,
                                                               .complementarityScore = 0.5,
                                                               .hybridizationEnergy = -18});

    EXPECT_TRUE(cluster1.overlaps(cluster2, 0));
}

TEST(InteractionClusterTest, OverlapsGrace) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{"record1", 0, GenomicStrand::FORWARD, 0, 10, 0.5, -18},
        RecordFragment{"record1", 1, GenomicStrand::FORWARD, 0, 10, 0.5, -18});
    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{"record2", 0, GenomicStrand::FORWARD, 11, 20, 0.5, -18},
        RecordFragment{"record2", 1, GenomicStrand::FORWARD, 0, 10, 0.5, -18});

    EXPECT_TRUE(cluster1.overlaps(cluster2, 1));
}

TEST(InteractionClusterTest, OverlapsNegativeGrace) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{"record1", 0, GenomicStrand::FORWARD, 0, 10, 0.5, -18},
        RecordFragment{"record1", 1, GenomicStrand::FORWARD, 0, 10, 0.5, -18});
    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{"record2", 0, GenomicStrand::FORWARD, 9, 20, 0.5, -18},
        RecordFragment{"record2", 1, GenomicStrand::FORWARD, 0, 10, 0.5, -18});

    EXPECT_TRUE(cluster1.overlaps(cluster2, -1));
}

TEST(InteractionClusterTest, NoOverlapsNegativeGrace) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{"record1", 0, GenomicStrand::FORWARD, 0, 10, 0.5, -18},
        RecordFragment{"record1", 1, GenomicStrand::FORWARD, 0, 10, 0.5, -18});
    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{"record2", 0, GenomicStrand::FORWARD, 9, 20, 0.5, -18},
        RecordFragment{"record2", 1, GenomicStrand::FORWARD, 0, 10, 0.5, -18});

    EXPECT_FALSE(cluster1.overlaps(cluster2, -2));
}

TEST(InteractionClusterTest, NoOverlapsFirstNotOverlapping) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{"record1", 0, GenomicStrand::FORWARD, 0, 10, 0.5, -18},
        RecordFragment{"record1", 1, GenomicStrand::FORWARD, 0, 10, 0.5, -18});
    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{"record2", 0, GenomicStrand::FORWARD, 11, 20, 0.5, -18},
        RecordFragment{"record2", 1, GenomicStrand::FORWARD, 0, 10, 0.5, -18});

    EXPECT_FALSE(cluster1.overlaps(cluster2, 0));
}

TEST(InteractionClusterTest, NoOverlapsSecondNotOverlapping) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{"record1", 0, GenomicStrand::FORWARD, 0, 10, 0.5, -18},
        RecordFragment{"record1", 1, GenomicStrand::FORWARD, 0, 10, 0.5, -18});
    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{"record2", 0, GenomicStrand::FORWARD, 0, 10, 0.5, -18},
        RecordFragment{"record2", 1, GenomicStrand::FORWARD, 11, 20, 0.5, -18});

    EXPECT_FALSE(cluster1.overlaps(cluster2, 0));
}

TEST(InteractionClusterTest, NoOverlapsDifferentReferenceID) {
    InteractionCluster cluster1 = InteractionCluster::fromRecordFragments(
        RecordFragment{"record1", 0, GenomicStrand::FORWARD, 0, 10, 0.5, -18},
        RecordFragment{"record1", 1, GenomicStrand::FORWARD, 0, 10, 0.5, -18});
    InteractionCluster cluster2 = InteractionCluster::fromRecordFragments(
        RecordFragment{"record2", 2, GenomicStrand::FORWARD, 0, 10, 0.5, -18},
        RecordFragment{"record2", 1, GenomicStrand::FORWARD, 0, 10, 0.5, -18});

    EXPECT_FALSE(cluster1.overlaps(cluster2, 0));
}

// InteractionCluster merge tests
TEST(InteractionClusterTest, BasicMerge) {
    InteractionCluster cluster1 =
        InteractionCluster::fromRecordFragments(RecordFragment{.recordID = "record1",
                                                               .referenceIDIndex = 0,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 5,
                                                               .end = 15,
                                                               .complementarityScore = 0.5,
                                                               .hybridizationEnergy = -18},
                                                RecordFragment{.recordID = "record1",
                                                               .referenceIDIndex = 0,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 20,
                                                               .end = 25,
                                                               .complementarityScore = 0.5,
                                                               .hybridizationEnergy = -18});
    InteractionCluster cluster2 =
        InteractionCluster::fromRecordFragments(RecordFragment{.recordID = "record2",
                                                               .referenceIDIndex = 0,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 0,
                                                               .end = 10,
                                                               .complementarityScore = 0.8,
                                                               .hybridizationEnergy = -15},
                                                RecordFragment{.recordID = "record2",
                                                               .referenceIDIndex = 0,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 20,
                                                               .end = 30,
                                                               .complementarityScore = 0.8,
                                                               .hybridizationEnergy = -15});

    cluster1.merge(cluster2);

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
    InteractionCluster cluster1 =
        InteractionCluster::fromRecordFragments(RecordFragment{.recordID = "record1",
                                                               .referenceIDIndex = 1,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 5,
                                                               .end = 15,
                                                               .complementarityScore = 0.5,
                                                               .hybridizationEnergy = -18},
                                                RecordFragment{.recordID = "record1",
                                                               .referenceIDIndex = 0,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 20,
                                                               .end = 25,
                                                               .complementarityScore = 0.5,
                                                               .hybridizationEnergy = -18});
    InteractionCluster cluster2 =
        InteractionCluster::fromRecordFragments(RecordFragment{.recordID = "record2",
                                                               .referenceIDIndex = 0,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 20,
                                                               .end = 30,
                                                               .complementarityScore = 0.8,
                                                               .hybridizationEnergy = -15},
                                                RecordFragment{.recordID = "record2",
                                                               .referenceIDIndex = 1,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 0,
                                                               .end = 10,
                                                               .complementarityScore = 0.8,
                                                               .hybridizationEnergy = -15});

    cluster1.merge(cluster2);

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
    InteractionCluster cluster1 =
        InteractionCluster::fromRecordFragments(RecordFragment{.recordID = "record1",
                                                               .referenceIDIndex = 0,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 5,
                                                               .end = 15,
                                                               .complementarityScore = 0.5,
                                                               .hybridizationEnergy = -18},
                                                RecordFragment{.recordID = "record1",
                                                               .referenceIDIndex = 0,
                                                               .strand = GenomicStrand::FORWARD,
                                                               .start = 20,
                                                               .end = 25,
                                                               .complementarityScore = 0.5,
                                                               .hybridizationEnergy = -18});
    InteractionCluster cluster2 =
        InteractionCluster::fromRecordFragments(RecordFragment{.recordID = "record2",
                                                               .referenceIDIndex = 0,
                                                               .strand = GenomicStrand::REVERSE,
                                                               .start = 0,
                                                               .end = 10,
                                                               .complementarityScore = 0.8,
                                                               .hybridizationEnergy = -15},
                                                RecordFragment{.recordID = "record2",
                                                               .referenceIDIndex = 0,
                                                               .strand = GenomicStrand::REVERSE,
                                                               .start = 20,
                                                               .end = 30,
                                                               .complementarityScore = 0.8,
                                                               .hybridizationEnergy = -15});

    EXPECT_DEATH(cluster1.merge(cluster2), ".*");
}
