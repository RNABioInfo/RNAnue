// NOLINTBEGIN

#pragma once

// Standard
#if __has_include(<execution>)
#include <execution>  // IWYU pragma: keep
#endif
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <unordered_set>
#include <vector>

// Internal
#include "VectorUtility.hpp"

/* Suppose there are N=2^(K+1)-1 sorted numbers in an array a[]. They
 * implicitly form a complete binary tree of height K+1. We consider leaves to
 * be at level 0. The binary tree has the following properties:
 *
 * 1. The lowest k-1 bits of nodes at level k are all 1. The k-th bit is 0.
 *    The first node at level k is indexed by 2^k-1. The root of the tree is
 *    indexed by 2^K-1.
 *
 * 2. For a node x at level k, its left child is x-2^(k-1) and the right child
 *    is x+2^(k-1).
 *
 * 3. For a node x at level k, it is a left child if its (k+1)-th bit is 0. Its
 *    parent node is x+2^k. Similarly, if the (k+1)-th bit is 1, x is a right
 *    child and its parent is x-2^k.
 *
 * 4. For a node x at level k, there are 2^(k+1)-1 nodes in the subtree
 *    descending from x, including x. The left-most leaf is x&~(2^k-1) (masking
 *    the lowest k bits to 0).
 *
 * When numbers can't fill a complete binary tree, the parent of a node may not
 * be present in the array. The implementation here still mimics a complete
 * tree, though getting the special casing right is a little complex. There may
 * be alternative solutions.
 *
 * As a sorted array can be considered as a binary search tree, we can
 * implement an interval tree on top of the idea. We only need to record, for
 * each node, the maximum value in the subtree descending from the node.
 */
template <typename S, typename T>  // "S" is a scalar type; "T" is the type of data associated with
                                   // each interval
class IITree {
   public:
    void add(const S& start, const S& end, const T& data) {
        assert(end > start && "Attempted adding invalid interval");
        intervalls.emplace_back(start, end, data);
    }

    /**
     * @brief Removes the interval at the specified index.
     *
     * This method erases the interval located at the given index from the internal
     * storage and reindexes the tree without sorting.
     *
     * @param index The index of the interval to remove.
     */
    void remove(size_t index) {
        intervalls.erase(intervalls.begin() + index);
        indexNoSort();
    }

    /**
     * @brief Removes multiple intervals specified by their indices.
     *
     * This method erases multiple intervals from the internal storage based on the
     * provided sorted vector of indices and reindexes the tree without sorting.
     *
     * @param indices A vector of indices of intervals to remove.
     */
    void remove(std::vector<size_t>& indices) {
        std::ranges::sort(indices);
        helper::erase_selected(intervalls, indices, intervalls.front());
        indexNoSort();
    }

    /**
     * @brief Removes multiple intervals specified by their indices.
     *
     * This method erases multiple intervals from the internal storage based on the
     * provided unordered_set of indices. It first converts the set to a sorted vector,
     * then removes the specified intervals and reindexes the tree without sorting.
     *
     * @param indices An unordered_set of indices of intervals to remove.
     */
    void remove(std::unordered_set<size_t>& indices) {
        std::vector<size_t> removalIndices{indices.begin(), indices.end()};
        std::ranges::sort(removalIndices);
        helper::erase_selected(intervalls, removalIndices, intervalls.front());
        indexNoSort();
    }

    void index() {
#ifdef __cpp_lib_execution
        std::sort(std::execution::par, intervalls.begin(), intervalls.end(), IntervalLess());
#else
        // Some compilers have not implemented execution policies, so don't.
        std::sort(intervalls.begin(), intervalls.end(), IntervalLess());
#endif

        max_level = index_core(intervalls);
    }

    void indexNoSort() { max_level = index_core(intervalls); }

    /**
     * @brief Checks for intervals that overlap with the specified [start, end) range (exclusive of
     * upper bound).
     *
     * This function traverses this interval tree to find intervals that overlap with
     * the input range [start, end). If overlaps are detected, their indices are
     * collected in the 'out' vector.
     *
     * @param start The start of the query interval.
     * @param end   The end of the query interval.
     * @param out   A reference to a vector that will receive the indices of all
     *              intervals that overlap the specified range.
     * @return True if any intervals overlap with the specified range; false otherwise.
     */
    auto overlap(const S& start, const S& end, std::vector<size_t>& out) const -> bool {
        assert(end > start);

        int t = 0;
        StackCell stack[64];
        out.clear();

        if (max_level < 0) {
            return false;
        }

        stack[t++] = StackCell(max_level, (1LL << max_level) - 1, 0);
        while (t) {
            StackCell z = stack[--t];
            if (z.k <= 3) {
                size_t i;
                size_t i0 = z.x >> z.k << z.k;
                size_t i1 = i0 + (1LL << (z.k + 1)) - 1;

                if (i1 >= intervalls.size()) {
                    i1 = intervalls.size();
                }

                for (i = i0; i < i1 && intervalls[i].start < end; ++i) {
                    if (start < intervalls[i].end) {
                        out.push_back(i);
                    }
                }
            } else if (z.w == 0) {
                size_t y = z.x - (1LL << (z.k - 1));
                stack[t++] = StackCell(z.k, z.x, 1);

                if (y >= intervalls.size() || intervalls[y].max > start) {
                    stack[t++] = StackCell(z.k - 1, y, 0);
                }
            } else if (z.x < intervalls.size() && intervalls[z.x].start < end) {
                if (start < intervalls[z.x].end) {
                    out.push_back(z.x);
                }

                stack[t++] = StackCell(z.k - 1, z.x + (1LL << (z.k - 1)), 0);
            }
        }

        return !out.empty();
    }

    /**
     * @brief Prints debug information about each interval in this data structure.
     */
    void printDebug() {
        for (size_t index = 0; index < intervalls.size(); ++index) {
            const auto& interval = intervalls[index];
            std::cout << "Index: " << index << "; " << interval.start << "-" << interval.end << "; "
                      << interval.data << std::endl;
        }
    }

    // Getters
    [[nodiscard]] constexpr auto size() const noexcept -> size_t { return intervalls.size(); }
    [[nodiscard]] constexpr auto getIntervalStart(size_t index) const noexcept -> const S& {
        return intervalls[index].start;
    }
    [[nodiscard]] constexpr auto getIntervalEnd(size_t index) const noexcept -> const S& {
        return intervalls[index].end;
    }
    [[nodiscard]] constexpr auto getData(size_t index) const noexcept -> const T& {
        return intervalls[index].data;
    }
    [[nodiscard]] constexpr auto getData(size_t index) noexcept -> T& {
        return intervalls[index].data;
    }

    // Setters
    void setIntervalStart(size_t index, const S& start) { intervalls[index].start = start; }
    void setIntervalEnd(size_t index, const S& end) { intervalls[index].end = end; }

   private:
    struct StackCell {
        size_t x{};
        int k{}, w{};
        StackCell() = default;
        StackCell(int k_, size_t x_, int w_) : x(x_), k(k_), w(w_) {}
    };

    struct Interval {
        S start, end, max;
        T data;
        Interval(const S& start, const S& end, const T& data)
            : start(start), end(end), max(end), data(data) {}

        Interval() = default;
    };

    struct IntervalLess {
        auto operator()(const Interval& lhs, const Interval& rhs) const -> bool {
            return lhs.start < rhs.start;
        }
    };

    std::vector<Interval> intervalls;
    int max_level{};

    auto index_core(std::vector<Interval>& intervalls) -> int {
        size_t i{};
        size_t last_i{};
        S last{};
        size_t k{};

        if (intervalls.empty()) {
            return -1;
        }

        for (i = 0; i < intervalls.size(); i += 2) {
            last_i = i, last = intervalls[i].max = intervalls[i].end;
        }

        for (k = 1; 1UL << k <= intervalls.size(); ++k) {
            size_t x = 1UL << (k - 1);
            size_t i0 = (x << 1) - 1;
            size_t step = x << 2;

            for (i = i0; i < intervalls.size(); i += step) {
                S el = intervalls[i - x].max;
                S er = i + x < intervalls.size() ? intervalls[i + x].max : last;
                S e = intervalls[i].end;
                e = e > el ? e : el;
                e = e > er ? e : er;
                intervalls[i].max = e;
            }

            last_i = last_i >> k & 1 ? last_i - x : last_i + x;

            if (last_i < intervalls.size() && intervalls[last_i].max > last) {
                last = intervalls[last_i].max;
            }
        }

        return k - 1;
    }

   public:
    [[nodiscard]] auto intervals() const -> const std::vector<Interval>& { return intervalls; }
};

// NOLINTEND
