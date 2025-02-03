#include <gtest/gtest.h>

#include <cstddef>
#include <string>
#include <utility>

// Internal
#include "IITree.hpp"
#include "annotation/GenomicFeatureTreeMerger.hpp"
#include "pipelines/dataTypes/GenomicFeature.hpp"
#include "pipelines/dataTypes/GenomicRegion.hpp"
#include "pipelines/dataTypes/GenomicStrand.hpp"
#include "pipelines/dataTypes/GenomicStrandSpecificity.hpp"

namespace {

using namespace annotation;  // NOLINT
using namespace dataTypes;   // NOLINT

class GenomicFeatureTreeMergerTest : public testing::Test {
   protected:
    // Helper method to create a GenomicFeature for convenience
    static auto makeFeature(int start, int end, GenomicStrand strand = GenomicStrand::FORWARD,
                            std::string featureType = "testFeature") -> GenomicFeature {
        // Set reference id to zero since IITree shpuld only contain only one reference ID
        GenomicRegion region(0, {.startPosition = start, .endPosition = end}, strand);
        return asGenomicFeature(region, std::move(featureType));
    }
};

// Test merging with an empty tree
TEST_F(GenomicFeatureTreeMergerTest, MergeEmptyTree) {
    IITree<int, GenomicFeature> featureTree;
    EXPECT_EQ(featureTree.size(), 0UL);

    FeatureMergingParameters params(GenomicStrandSpecificity::SPECIFIC, /*mergingTolerance=*/10);
    GenomicFeatureTreeMerger merger(params);
    merger.merge(featureTree);

    // Tree should remain empty
    EXPECT_EQ(featureTree.size(), 0UL);
}

// Test merging a single feature: nothing to merge
TEST_F(GenomicFeatureTreeMergerTest, MergeSingleFeature) {
    IITree<int, GenomicFeature> featureTree;
    featureTree.add(1, 5, makeFeature(1, 5));
    featureTree.index();
    EXPECT_EQ(featureTree.size(), 1UL);

    FeatureMergingParameters params(GenomicStrandSpecificity::SPECIFIC, /*mergingTolerance=*/10);
    GenomicFeatureTreeMerger merger(params);

    merger.merge(featureTree);
    // Should remain a single feature
    EXPECT_EQ(featureTree.size(), 1UL);
    EXPECT_EQ(featureTree.getIntervalStart(0), 1);
    EXPECT_EQ(featureTree.getIntervalEnd(0), 5);
}

// Test merging with no overlapping intervals
TEST_F(GenomicFeatureTreeMergerTest, MergeNonOverlapping) {
    IITree<int, GenomicFeature> featureTree;
    featureTree.add(0, 5, makeFeature(0, 5));
    featureTree.add(10, 15, makeFeature(10, 15));
    featureTree.add(20, 25, makeFeature(20, 25));
    featureTree.index();
    EXPECT_EQ(featureTree.size(), 3UL);

    FeatureMergingParameters params(GenomicStrandSpecificity::SPECIFIC, /*mergingTolerance=*/0);
    GenomicFeatureTreeMerger merger(params);
    merger.merge(featureTree);

    // None of the features should merge because there's no overlap (tolerance=0).
    EXPECT_EQ(featureTree.size(), 3UL);

    // Ensure intervals are still as originally set.
    EXPECT_EQ(featureTree.getIntervalStart(0), 0);
    EXPECT_EQ(featureTree.getIntervalEnd(0), 5);
    EXPECT_EQ(featureTree.getIntervalStart(1), 10);
    EXPECT_EQ(featureTree.getIntervalEnd(1), 15);
    EXPECT_EQ(featureTree.getIntervalStart(2), 20);
    EXPECT_EQ(featureTree.getIntervalEnd(2), 25);
}

// Test merging fully overlapping intervals with zero tolerance
TEST_F(GenomicFeatureTreeMergerTest, MergeFullyOverlappingZeroTolerance) {
    IITree<int, GenomicFeature> featureTree;
    // Intervals [1, 10), [2, 8), [3, 9) all overlap each other
    featureTree.add(1, 10, makeFeature(1, 10));
    featureTree.add(2, 8, makeFeature(2, 8));
    featureTree.add(3, 9, makeFeature(3, 9));
    // Non-overlapped interval
    featureTree.add(20, 30, makeFeature(20, 30));

    featureTree.index();
    EXPECT_EQ(featureTree.size(), 4UL);

    FeatureMergingParameters params(GenomicStrandSpecificity::SPECIFIC, /*mergingTolerance=*/0);
    GenomicFeatureTreeMerger merger(params);
    merger.merge(featureTree);

    // Expect the first three intervals to merge into one, plus the non-overlapping.
    EXPECT_EQ(featureTree.size(), 2UL);

    // The merged interval should reflect the bounding region of [1, 10).
    EXPECT_EQ(featureTree.getIntervalStart(0), 1);
    EXPECT_EQ(featureTree.getIntervalEnd(0), 10);

    // The second interval remains the non-overlapped one.
    EXPECT_EQ(featureTree.getIntervalStart(1), 20);
    EXPECT_EQ(featureTree.getIntervalEnd(1), 30);
}

// Test merging partially overlapping intervals with zero tolerance
TEST_F(GenomicFeatureTreeMergerTest, MergePartialOverlapZeroTolerance) {
    IITree<int, GenomicFeature> featureTree;
    // Intervals [0,5) and [5,10) only "touch" at boundary 5 => with zero tolerance they do NOT
    // merge.
    featureTree.add(0, 5, makeFeature(0, 5));
    featureTree.add(5, 10, makeFeature(5, 10));

    featureTree.index();

    EXPECT_EQ(featureTree.size(), 2UL);

    FeatureMergingParameters params(GenomicStrandSpecificity::SPECIFIC, /*mergingTolerance=*/0);
    GenomicFeatureTreeMerger merger(params);
    merger.merge(featureTree);

    // They remain two intervals because there's no overlapping interior.
    EXPECT_EQ(featureTree.size(), 2UL);
    EXPECT_EQ(featureTree.getIntervalStart(0), 0);
    EXPECT_EQ(featureTree.getIntervalEnd(0), 5);
    EXPECT_EQ(featureTree.getIntervalStart(1), 5);
    EXPECT_EQ(featureTree.getIntervalEnd(1), 10);
}

// Test merging intervals with a nonzero tolerance
TEST_F(GenomicFeatureTreeMergerTest, MergeWithTolerance) {
    IITree<int, GenomicFeature> featureTree;
    // Intervals [0,5) and [6,10) have a gap of 1. Tolerance=2 => they should merge.
    featureTree.add(0, 5, makeFeature(0, 5));
    featureTree.add(6, 10, makeFeature(6, 10));

    featureTree.index();

    EXPECT_EQ(featureTree.size(), 2UL);

    FeatureMergingParameters params(GenomicStrandSpecificity::SPECIFIC, /*mergingTolerance=*/2);
    GenomicFeatureTreeMerger merger(params);
    merger.merge(featureTree);

    // Expect these two intervals to merge into a single interval [0,10)
    EXPECT_EQ(featureTree.size(), 1UL);
    EXPECT_EQ(featureTree.getIntervalStart(0), 0);
    EXPECT_EQ(featureTree.getIntervalEnd(0), 10);
}

// Test merging intervals in UNSPECIFIC strand mode
TEST_F(GenomicFeatureTreeMergerTest, MergeUnspecificStrand) {
    IITree<int, GenomicFeature> featureTree;
    // Forward feature
    featureTree.add(1, 5, makeFeature(1, 5, GenomicStrand::FORWARD));
    // Reverse feature partially overlapping
    featureTree.add(4, 9, makeFeature(4, 9, GenomicStrand::REVERSE));
    // Reverse feature not overlapping
    featureTree.add(20, 25, makeFeature(20, 25, GenomicStrand::REVERSE));

    featureTree.index();

    EXPECT_EQ(featureTree.size(), 3UL);

    // In UNSPECIFIC mode, merges happen regardless of strand differences
    FeatureMergingParameters params(GenomicStrandSpecificity::UNSPECIFIC, /*mergingTolerance=*/0);
    GenomicFeatureTreeMerger merger(params);
    merger.merge(featureTree);

    // The first two intervals overlap in range: [1,5), [4,9) => merges into [1,9), ignoring strand
    // The non-overlapping [20,25) remains
    EXPECT_EQ(featureTree.size(), 2UL);
    EXPECT_EQ(featureTree.getIntervalStart(0), 1);
    EXPECT_EQ(featureTree.getIntervalEnd(0), 9);
    EXPECT_EQ(featureTree.getIntervalStart(1), 20);
    EXPECT_EQ(featureTree.getIntervalEnd(1), 25);

    // The merged interval should have .strand = NONE
    // Because the code sets it to NONE in UNSPECIFIC merges:
    EXPECT_EQ(featureTree.getData(0).genomicRegion.getStrand(), GenomicStrand::NONE);
}

// Test merging intervals in SPECIFIC strand mode
TEST_F(GenomicFeatureTreeMergerTest, MergeSpecificStrand) {
    IITree<int, GenomicFeature> featureTree;
    // Forward feature
    featureTree.add(1, 5, makeFeature(1, 5, GenomicStrand::FORWARD));
    // Reverse feature overlaps by position but not by strand
    featureTree.add(2, 6, makeFeature(2, 6, GenomicStrand::REVERSE));
    // Additional forward feature overlapping in range with the first
    featureTree.add(4, 7, makeFeature(4, 7, GenomicStrand::FORWARD));

    featureTree.index();

    EXPECT_EQ(featureTree.size(), 3UL);

    // SPECIFIC mode => only merges if orientation has the same strand
    FeatureMergingParameters params(GenomicStrandSpecificity::SPECIFIC, /*mergingTolerance=*/0);
    GenomicFeatureTreeMerger merger(params);
    merger.merge(featureTree);

    // Expect the two forward features [1,5) and [4,7) => merges into [1,7)
    // The reverse feature remains separate
    EXPECT_EQ(featureTree.size(), 2UL);

    // Check that the merged forward interval is [1,7)
    // The next interval is the reverse one [2,6)
    // Because we always keep the "baseFeature" in place, the first item in the tree will be the
    // merged forward region
    EXPECT_EQ(featureTree.getIntervalStart(0), 1);
    EXPECT_EQ(featureTree.getIntervalEnd(0), 7);
    EXPECT_EQ(featureTree.getIntervalStart(1), 2);
    EXPECT_EQ(featureTree.getIntervalEnd(1), 6);

    // Check strands
    EXPECT_EQ(featureTree.getData(0).genomicRegion.getStrand(), GenomicStrand::FORWARD);
    EXPECT_EQ(featureTree.getData(1).genomicRegion.getStrand(), GenomicStrand::REVERSE);
}

// Test that multiple merges happen iteratively
TEST_F(GenomicFeatureTreeMergerTest, MergeMultipleStages) {
    IITree<int, GenomicFeature> featureTree;
    featureTree.add(1, 3, makeFeature(1, 3));    // Overlaps with next
    featureTree.add(2, 6, makeFeature(2, 6));    // Overlaps with above & next
    featureTree.add(5, 8, makeFeature(5, 8));    // Overlaps with above & next
    featureTree.add(7, 10, makeFeature(7, 10));  // Overlaps with above
    // Another disjoint one
    featureTree.add(20, 25, makeFeature(20, 25));

    featureTree.index();

    EXPECT_EQ(featureTree.size(), 5UL);

    FeatureMergingParameters params(GenomicStrandSpecificity::SPECIFIC, /*mergingTolerance=*/0);
    GenomicFeatureTreeMerger merger(params);
    merger.merge(featureTree);

    // The first 4 features all chain-merge into one => [1,10)
    // The last remains separate => [20,25)
    EXPECT_EQ(featureTree.size(), 2UL);
    EXPECT_EQ(featureTree.getIntervalStart(0), 1);
    EXPECT_EQ(featureTree.getIntervalEnd(0), 10);
    EXPECT_EQ(featureTree.getIntervalStart(1), 20);
    EXPECT_EQ(featureTree.getIntervalEnd(1), 25);
}

// Test that removal indices are correctly handled inside merge
TEST_F(GenomicFeatureTreeMergerTest, MergeAndCheckRemovedIndices) {
    IITree<int, GenomicFeature> featureTree;
    // Intervals that do NOT all overlap among themselves, but partial
    // This ensures some intervals get removed while others remain.
    // Overlapping chain: [1,4), [2,3)
    featureTree.add(1, 4, makeFeature(1, 4));
    featureTree.add(2, 3, makeFeature(2, 3));

    // Next chain: [7,8), [6,9), [8,9)
    featureTree.add(7, 8, makeFeature(7, 8));
    featureTree.add(6, 9, makeFeature(6, 9));
    featureTree.add(8, 9, makeFeature(8, 9));

    featureTree.index();

    EXPECT_EQ(featureTree.size(), 5UL);

    FeatureMergingParameters params(GenomicStrandSpecificity::SPECIFIC, /*mergingTolerance=*/0);
    GenomicFeatureTreeMerger merger(params);
    merger.merge(featureTree);

    // [1,4), [2,3) => merge => [1,4)
    // [7,8), [6,9), [8,9) => merge => [6,9)
    // => final 2 intervals: [1,4), [6,9)
    EXPECT_EQ(featureTree.size(), 2UL);
    EXPECT_EQ(featureTree.getIntervalStart(0), 1);
    EXPECT_EQ(featureTree.getIntervalEnd(0), 4);
    EXPECT_EQ(featureTree.getIntervalStart(1), 6);
    EXPECT_EQ(featureTree.getIntervalEnd(1), 9);
}

}  // end anonymous namespace
