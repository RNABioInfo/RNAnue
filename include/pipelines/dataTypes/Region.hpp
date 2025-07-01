#pragma once

// Standard
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace dataTypes {

/**
 * @brief Represents region.
 *
 * Regions are always 0-based, half-open intervals.
 */
struct Region {
    int32_t startPosition;
    int32_t endPosition;

    [[nodiscard]] constexpr auto length() const -> int32_t { return endPosition - startPosition; }

    /**
     * @brief Computes the overlap between this region and another region.
     *
     * This function calculates the number of bases in which the two regions overlap.
     * If there is no overlap, it returns zero.
     *
     * @param other The other region to compute the overlap with.
     * @return The size of the overlapping interval between the two regions.
     */
    [[nodiscard]] constexpr auto overlap(const Region& other) const noexcept -> size_t {
        return std::max(0, std::min(endPosition, other.endPosition) -
                               std::max(startPosition, other.startPosition));
    }

    /**
     * @brief Checks whether two regions overlap.
     * @param other The other region.
     */
    [[nodiscard]] auto overlaps(const Region& other) const noexcept -> bool {
        return endPosition > other.startPosition && other.endPosition > startPosition;
    }

    /**
     * @brief Checks whether two regions overlap.
     *
     * Overlap means at least one nucleotide overlaps e.g. blunt ended regions with tolerance 0 do
     * not overlap.
     *
     * @param other The other region.
     * @param tolerance The tolerance to apply to the overlap calculation. 0 means at least one base
     * overlap, negative requires more bases to overlap and positive allows distance.
     */
    [[nodiscard]] auto overlapsWithTolerance(const Region& other,
                                             const int tolerance) const noexcept -> bool {
        return (startPosition <= other.startPosition)
                   ? (other.startPosition - endPosition < tolerance)
                   : (startPosition - other.endPosition < tolerance);
    }

    /**
     * @brief Checks if the overlap between this region and another region meets or exceeds a
     * specified fraction of the smaller region's length.
     *
     * This method computes the overlap between this region and another region, then compares it
     * against the length of the shorter region. If the overlap is at least the specified fraction
     * of that shorter region, the method returns true.
     *
     * @param other The other region to compare against.
     * @param overlapFraction The fraction of the smaller region's length that must be overlapped.
     * @return True if the overlap meets or exceeds the specified fraction of the smaller region,
     * otherwise false.
     */
    [[nodiscard]] auto overlapsWithShortestRegionFraction(
        const Region& other, const float overlapFraction) const noexcept -> bool {
        assert(FP_ZERO != std::fpclassify(overlapFraction) &&
               "Overlap fraction must be larger then 0.0");

        const size_t shortestLength = std::min(length(), other.length());
        if (shortestLength == 0) {
            return false;
        }

        const auto overlapLength = this->overlap(other);
        if (overlapLength == 0) {
            return false;
        }

        const float fraction =
            static_cast<float>(overlapLength) / static_cast<float>(shortestLength);
        return fraction >= overlapFraction;
    }

    /**
     * @brief Merges this region with another region.
     *
     * Sets the startPosition to the minimum of the two regions' start values,
     * and the endPosition to the maximum of the two regions' end values.
     *
     * @param other The other region to merge with.
     */
    auto merge(const Region& other) {
        startPosition = std::min(startPosition, other.startPosition);
        endPosition = std::max(endPosition, other.endPosition);
    }

    /**
     * @brief Expands the region by a given number of bases.
     * Region will not be be expanded below 0.
     * @param amount The number of bases to expand the region by.
     * @return Reference to this object.
     */
    auto expand(const size_t amount) -> Region& {
        startPosition = std::max<int32_t>(0, startPosition - static_cast<int32_t>(amount));
        endPosition += static_cast<int32_t>(amount);
        return *this;
    };

    /**
     * @brief Checks if the given position is contained within the region.
     *
     * Determines whether the specified position lies between the region's
     * startPosition (inclusive) and endPosition (exclusive).
     *
     * @param position The position to be checked.
     * @return True if the position is contained, otherwise false.
     */
    [[nodiscard]] auto contains(const int32_t position) const -> bool {
        return startPosition <= position && position <= endPosition;
    }

    // Operators
    [[nodiscard]] auto operator==(const Region& other) const -> bool {
        return startPosition == other.startPosition && endPosition == other.endPosition;
    }
};

}  // namespace dataTypes
