// NOLINTBEGIN
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <string>
#include <unordered_set>
#include <vector>

#include "gtest/gtest.h"

// Internal
#include "IITree.hpp"

class IITreeTest : public testing::Test {
   protected:
    IITree<int, std::string> tree;
};

// Test adding intervals and checking size
TEST_F(IITreeTest, AddIntervals) {
    EXPECT_EQ(tree.size(), 0UL);

    tree.add(1, 5, "Interval1");
    EXPECT_EQ(tree.size(), 1UL);
    EXPECT_EQ(tree.getIntervalStart(0), 1);
    EXPECT_EQ(tree.getIntervalEnd(0), 5);
    EXPECT_EQ(tree.getData(0), "Interval1");

    tree.add(3, 7, "Interval2");
    EXPECT_EQ(tree.size(), 2UL);
    EXPECT_EQ(tree.getIntervalStart(1), 3);
    EXPECT_EQ(tree.getIntervalEnd(1), 7);
    EXPECT_EQ(tree.getData(1), "Interval2");
}

// Test removing a single interval by index
TEST_F(IITreeTest, RemoveSingleInterval) {
    tree.add(1, 5, "Interval1");
    tree.add(3, 7, "Interval2");
    tree.add(6, 10, "Interval3");
    EXPECT_EQ(tree.size(), 3UL);

    tree.remove(1);  // Remove "Interval2"
    EXPECT_EQ(tree.size(), 2UL);
    EXPECT_EQ(tree.getIntervalStart(0), 1);
    EXPECT_EQ(tree.getIntervalEnd(0), 5);
    EXPECT_EQ(tree.getData(0), "Interval1");
    EXPECT_EQ(tree.getIntervalStart(1), 6);
    EXPECT_EQ(tree.getIntervalEnd(1), 10);
    EXPECT_EQ(tree.getData(1), "Interval3");
}

// Test removing multiple intervals using a vector of indices
TEST_F(IITreeTest, RemoveMultipleIntervals_Vector) {
    tree.add(1, 5, "Interval1");
    tree.add(3, 7, "Interval2");
    tree.add(6, 10, "Interval3");
    tree.add(8, 12, "Interval4");
    EXPECT_EQ(tree.size(), 4u);

    std::vector<size_t> indices = {1, 3};  // Remove "Interval2" and "Interval4"
    tree.remove(indices);
    EXPECT_EQ(tree.size(), 2u);
    EXPECT_EQ(tree.getIntervalStart(0), 1);
    EXPECT_EQ(tree.getIntervalEnd(0), 5);
    EXPECT_EQ(tree.getData(0), "Interval1");
    EXPECT_EQ(tree.getIntervalStart(1), 6);
    EXPECT_EQ(tree.getIntervalEnd(1), 10);
    EXPECT_EQ(tree.getData(1), "Interval3");
}

// Test removing multiple intervals using an unordered_set of indices
TEST_F(IITreeTest, RemoveMultipleIntervals_UnorderedSet) {
    tree.add(1, 5, "Interval1");
    tree.add(3, 7, "Interval2");
    tree.add(6, 10, "Interval3");
    tree.add(8, 12, "Interval4");
    EXPECT_EQ(tree.size(), 4UL);

    std::unordered_set<size_t> indices = {0, 2};  // Remove "Interval1" and "Interval3"
    tree.remove(indices);
    EXPECT_EQ(tree.size(), 2UL);
    EXPECT_EQ(tree.getIntervalStart(0), 3);
    EXPECT_EQ(tree.getIntervalEnd(0), 7);
    EXPECT_EQ(tree.getData(0), "Interval2");
    EXPECT_EQ(tree.getIntervalStart(1), 8);
    EXPECT_EQ(tree.getIntervalEnd(1), 12);
    EXPECT_EQ(tree.getData(1), "Interval4");
}

// Test indexing the tree and verifying functionality
TEST_F(IITreeTest, IndexTree) {
    // Add intervals in random order
    tree.add(10, 20, "Interval1");
    tree.add(5, 15, "Interval2");
    tree.add(20, 30, "Interval3");
    tree.add(25, 35, "Interval4");
    tree.add(0, 4, "Interval5");
    EXPECT_EQ(tree.size(), 5UL);

    // Index the tree
    tree.index();
    // Since index_core returns max_level, but it's private, we can't directly check it.
    // Instead, perform an overlap query which relies on correct indexing.

    std::vector<size_t> overlaps;
    bool found = tree.overlap(3, 6, overlaps);
    EXPECT_TRUE(found);
    // Expected overlapping intervals: Interval2 ([5,15)), Interval5 ([0,4)) overlaps [3,6)
    // However, Interval1 starts at 10 which is > 6, so only Interval2 and Interval5
    EXPECT_EQ(overlaps.size(), 2UL);
    EXPECT_TRUE((overlaps == std::vector<size_t>{0, 1}));
}

// Test indexing without sorting
TEST_F(IITreeTest, IndexNoSort) {
    // Add intervals already sorted
    tree.add(0, 5, "Interval1");
    tree.add(6, 10, "Interval2");
    tree.add(11, 15, "Interval3");
    EXPECT_EQ(tree.size(), 3UL);

    tree.indexNoSort();

    std::vector<size_t> overlaps;
    bool found = tree.overlap(7, 12, overlaps);
    EXPECT_TRUE(found);
    // Overlapping intervals: Interval2 and Interval3
    EXPECT_EQ(overlaps.size(), 2UL);
    EXPECT_TRUE((overlaps == std::vector<size_t>{1, 2}));
}

// Test overlap queries
TEST_F(IITreeTest, OverlapQueries) {
    // Add multiple intervals
    tree.add(1, 5, "Interval1");
    tree.add(3, 7, "Interval2");
    tree.add(6, 10, "Interval3");
    tree.add(8, 12, "Interval4");
    tree.add(11, 15, "Interval5");
    tree.index();

    // Overlap with [4,9)
    std::vector<size_t> overlaps;
    bool found = tree.overlap(4, 9, overlaps);
    EXPECT_TRUE(found);
    // Overlapping intervals: Interval1 ([1,5)), Interval2 ([3,7)), Interval3 ([6,10)), Interval4
    // ([8,12))
    EXPECT_EQ(overlaps.size(), 4UL);
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 0) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 1) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 2) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 3) != overlaps.end()));

    // Overlap with [0,1)
    overlaps.clear();
    found = tree.overlap(0, 1, overlaps);
    EXPECT_FALSE(found);
    EXPECT_TRUE(overlaps.empty());

    // Overlap with [10, 14)
    overlaps.clear();
    found = tree.overlap(10, 14, overlaps);
    EXPECT_TRUE(found);
    // Overlapping intervals: Interval4 ([8,12)), Interval5 ([11,15))
    EXPECT_EQ(overlaps.size(), 2UL);
    EXPECT_TRUE((overlaps == std::vector<size_t>{3, 4}));

    // Overlap with [15, 20)
    overlaps.clear();
    found = tree.overlap(15, 20, overlaps);
    EXPECT_FALSE(found);
    EXPECT_TRUE(overlaps.empty());
}

// Test getters and setters
TEST_F(IITreeTest, GettersSetters) {
    tree.add(2, 6, "Interval1");
    tree.add(5, 9, "Interval2");
    EXPECT_EQ(tree.size(), 2UL);

    // Test getters
    EXPECT_EQ(tree.getIntervalStart(0), 2);
    EXPECT_EQ(tree.getIntervalEnd(0), 6);
    EXPECT_EQ(tree.getData(0), "Interval1");
    EXPECT_EQ(tree.getIntervalStart(1), 5);
    EXPECT_EQ(tree.getIntervalEnd(1), 9);
    EXPECT_EQ(tree.getData(1), "Interval2");

    // Test setters
    tree.setIntervalStart(0, 3);
    tree.setIntervalEnd(0, 7);
    tree.getData(0) = "Interval1_Modified";

    EXPECT_EQ(tree.getIntervalStart(0), 3);
    EXPECT_EQ(tree.getIntervalEnd(0), 7);
    EXPECT_EQ(tree.getData(0), "Interval1_Modified");
}

// Test removing intervals and ensuring indices are updated correctly
TEST_F(IITreeTest, RemoveAndValidateIndices) {
    // Add intervals
    tree.add(1, 4, "A");
    tree.add(2, 5, "B");
    tree.add(3, 6, "C");
    tree.add(4, 7, "D");
    tree.add(5, 8, "E");
    EXPECT_EQ(tree.size(), 5UL);

    // Remove middle interval
    tree.remove(2);  // Remove "C"
    EXPECT_EQ(tree.size(), 4UL);
    EXPECT_EQ(tree.getIntervalStart(0), 1);
    EXPECT_EQ(tree.getIntervalEnd(0), 4);
    EXPECT_EQ(tree.getData(0), "A");
    EXPECT_EQ(tree.getIntervalStart(1), 2);
    EXPECT_EQ(tree.getIntervalEnd(1), 5);
    EXPECT_EQ(tree.getData(1), "B");
    EXPECT_EQ(tree.getIntervalStart(2), 4);
    EXPECT_EQ(tree.getIntervalEnd(2), 7);
    EXPECT_EQ(tree.getData(2), "D");
    EXPECT_EQ(tree.getIntervalStart(3), 5);
    EXPECT_EQ(tree.getIntervalEnd(3), 8);
    EXPECT_EQ(tree.getData(3), "E");

    // Remove first and last
    std::vector<size_t> indices = {0, 3};
    tree.remove(indices);
    EXPECT_EQ(tree.size(), 2UL);
    EXPECT_EQ(tree.getIntervalStart(0), 2);
    EXPECT_EQ(tree.getIntervalEnd(0), 5);
    EXPECT_EQ(tree.getData(0), "B");
    EXPECT_EQ(tree.getIntervalStart(1), 4);
    EXPECT_EQ(tree.getIntervalEnd(1), 7);
    EXPECT_EQ(tree.getData(1), "D");
}

// Test overlap queries after modifications
TEST_F(IITreeTest, OverlapAfterModifications) {
    // Add and index intervals
    tree.add(10, 20, "Interval1");
    tree.add(15, 25, "Interval2");
    tree.add(20, 30, "Interval3");
    tree.index();

    // Initial overlap
    std::vector<size_t> overlaps;
    bool found = tree.overlap(18, 22, overlaps);
    EXPECT_TRUE(found);
    EXPECT_EQ(overlaps.size(), 3UL);
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 0) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 1) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 2) != overlaps.end()));

    // Remove one interval and re-query
    tree.remove(1);  // Remove "Interval2"
    overlaps.clear();
    found = tree.overlap(18, 22, overlaps);
    EXPECT_TRUE(found);
    EXPECT_EQ(overlaps.size(), 2UL);
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 0) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 1) != overlaps.end()));
}

// Test handling of empty tree
TEST_F(IITreeTest, EmptyTree) {
    EXPECT_EQ(tree.size(), 0UL);

    std::vector<size_t> overlaps;
    bool found = tree.overlap(0, 10, overlaps);
    EXPECT_FALSE(found);
    EXPECT_TRUE(overlaps.empty());
}

// Test overlapping at boundaries
TEST_F(IITreeTest, OverlapBoundaries) {
    tree.add(5, 10, "Interval1");
    tree.add(10, 15, "Interval2");
    tree.add(15, 20, "Interval3");
    tree.index();

    std::vector<size_t> overlaps;
    bool found;

    // Query that touches the end of Interval1 and start of Interval2
    overlaps.clear();
    EXPECT_DEATH(tree.overlap(10, 10, overlaps), "");
    EXPECT_TRUE(overlaps.empty());

    // Query that includes the end of Interval1
    overlaps.clear();
    found = tree.overlap(9, 10, overlaps);
    EXPECT_TRUE(found);
    EXPECT_EQ(overlaps.size(), 1UL);
    EXPECT_EQ(overlaps[0], 0UL);

    // Query that includes the start of Interval2
    overlaps.clear();
    found = tree.overlap(10, 11, overlaps);
    EXPECT_TRUE(found);
    EXPECT_EQ(overlaps.size(), 1UL);
    EXPECT_EQ(overlaps[0], 1UL);
}

// Test multiple overlapping intervals
TEST_F(IITreeTest, MultipleOverlaps) {
    // Add overlapping intervals
    tree.add(1, 10, "A");
    tree.add(2, 9, "B");
    tree.add(3, 8, "C");
    tree.add(4, 7, "D");
    tree.add(5, 6, "E");
    tree.index();

    // Query overlapping with [4,5)
    std::vector<size_t> overlaps;
    bool found = tree.overlap(4, 5, overlaps);
    EXPECT_TRUE(found);
    EXPECT_EQ(overlaps.size(), 4UL);
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 0) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 1) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 2) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 3) != overlaps.end()));
}

// Test large number of intervals for performance and correctness
TEST_F(IITreeTest, LargeNumberOfIntervals) {
    const int num_intervals = 1000;
    for (int i = 0; i < num_intervals; ++i) {
        tree.add(i, i + 10, "Interval" + std::to_string(i));
    }
    tree.index();

    // Perform multiple overlap queries
    for (int q = 0; q < 100; ++q) {
        int start = q * 10;
        int end = start + 5;
        std::vector<size_t> overlaps;
        bool found = tree.overlap(start, end, overlaps);
        EXPECT_TRUE(found);
        // Ensure that all overlapped intervals indeed overlap the query
        for (size_t idx : overlaps) {
            EXPECT_TRUE(tree.getIntervalStart(idx) < end);
            EXPECT_TRUE(start < tree.getIntervalEnd(idx));
        }
    }

    // Check size
    EXPECT_EQ(tree.size(), static_cast<size_t>(num_intervals));
}

// Test overlapping all intervals
TEST_F(IITreeTest, OverlapAll) {
    // All intervals overlap with [0, 100)
    for (int i = 0; i < 50; ++i) {
        tree.add(i, 100, "Interval" + std::to_string(i));
    }
    tree.index();

    std::vector<size_t> overlaps;
    bool found = tree.overlap(50, 60, overlaps);
    EXPECT_TRUE(found);
    EXPECT_EQ(overlaps.size(), 50UL);
}

// Test non-overlapping consecutive intervals
TEST_F(IITreeTest, NonOverlappingConsecutive) {
    // Create consecutive non-overlapping intervals
    for (int i = 0; i < 100; ++i) {
        tree.add(i * 10, (i + 1) * 10, "Interval" + std::to_string(i));
    }
    tree.index();

    // Query each interval
    for (int i = 0; i < 100; ++i) {
        std::vector<size_t> overlaps;
        bool found = tree.overlap(i * 10, (i + 1) * 10, overlaps);
        EXPECT_TRUE(found);
        EXPECT_EQ(overlaps.size(), 1UL);
        EXPECT_EQ(overlaps[0], static_cast<size_t>(i));
    }

    // Query overlapping two intervals
    std::vector<size_t> overlaps;
    bool found = tree.overlap(15, 25, overlaps);
    EXPECT_TRUE(found);
    EXPECT_EQ(overlaps.size(), 2UL);
    EXPECT_EQ(overlaps[0], 1UL);
    EXPECT_EQ(overlaps[1], 2UL);
}

// Test overlapping with all intervals start before the query
TEST_F(IITreeTest, OverlapAllStartBefore) {
    // All intervals start before 50 but overlap into it
    for (int i = 0; i < 100; ++i) {
        tree.add(i, i + 100, "Interval" + std::to_string(i));
    }
    tree.index();

    std::vector<size_t> overlaps;
    bool found = tree.overlap(50, 150, overlaps);
    EXPECT_TRUE(found);
    EXPECT_EQ(overlaps.size(), 100UL);
    for (size_t idx : overlaps) {
        EXPECT_TRUE(tree.getIntervalStart(idx) < 150);
        EXPECT_TRUE(50 < tree.getIntervalEnd(idx));
    }
}

// Test intervals with the same start and different ends
TEST_F(IITreeTest, SameStartDifferentEnds) {
    tree.add(10, 20, "A");
    tree.add(10, 25, "B");
    tree.add(10, 15, "C");
    tree.add(10, 30, "D");
    tree.index();

    std::vector<size_t> overlaps;
    bool found = tree.overlap(10, 20, overlaps);
    EXPECT_TRUE(found);
    EXPECT_EQ(overlaps.size(), 4UL);
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 0) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 1) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 2) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 3) != overlaps.end()));
}

// Test intervals with same end and different starts
TEST_F(IITreeTest, SameEndDifferentStarts) {
    tree.add(5, 15, "A");
    tree.add(10, 15, "B");
    tree.add(12, 15, "C");
    tree.add(0, 15, "D");
    tree.index();

    std::vector<size_t> overlaps;
    bool found = tree.overlap(10, 15, overlaps);
    EXPECT_TRUE(found);
    EXPECT_EQ(overlaps.size(), 4UL);
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 0) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 1) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 2) != overlaps.end()));
    EXPECT_TRUE((std::find(overlaps.begin(), overlaps.end(), 3) != overlaps.end()));
}

// NOLINTEND
