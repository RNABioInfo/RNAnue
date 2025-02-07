#pragma once

// Standard
#include <zconf.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <optional>
#include <ostream>

// Internal
#include "GenomicOrientation.hpp"
#include "GenomicStrand.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "Region.hpp"
#include "SamRecord.hpp"

namespace dataTypes {
/**
 * @brief Represents a genomic region.
 *
 * Genomic regions are always 0-based, half-open intervals.
 */
struct GenomicRegion {
    /**
     * @brief Constructs a GenomicRegion object.
     * @param referenceID The reference ID of the genomic region.
     * @param region The region of the genomic region.
     * @param strand The strand of the genomic region (optional).
     */
    constexpr GenomicRegion(int referenceID, Region region,
                            GenomicStrand strand = GenomicStrand::NONE)
        : referenceIDIndex(referenceID), region(region), strand(strand) {};

    // Getters
    [[nodiscard]] constexpr auto getReferenceIDIndex() const { return referenceIDIndex; }
    [[nodiscard]] constexpr auto getStrand() const { return strand; }
    [[nodiscard]] constexpr auto getStart() const { return region.startPosition; }
    [[nodiscard]] constexpr auto getEnd() const { return region.endPosition; }

    // Setters
    void setReferenceIDIndex(int referenceID) noexcept { referenceIDIndex = referenceID; }
    void setStrand(GenomicStrand strand) noexcept { this->strand = strand; }
    void setStart(int start) noexcept { region.startPosition = start; }
    void setEnd(int end) noexcept { region.endPosition = end; }

    /**
     * @brief Creates a GenomicRegion object from a SamRecord.
     * @param record The SamRecord object.
     * @param referenceIDs The deque of reference IDs.
     * @return An optional GenomicRegion object.
     */
    static auto fromSamRecord(const SamRecord& record) -> std::optional<GenomicRegion>;

    /**
     * @brief Returns the length of this genomic region.
     *
     * Calculates and returns the difference between the end
     * and the start positions of the underlying region.
     *
     * @return The size of the region in base pairs.
     */
    [[nodiscard]] constexpr auto length() const noexcept -> size_t { return getEnd() - getStart(); }

    /**
     * @brief Expands the region into both direction by the specified amount
     * Region will not be expanded beyond the start of the reference
     * @param amount The amount to expand the region
     * @return Reference to the GenomicRegion object
     */
    auto expand(const size_t amount) -> GenomicRegion& {
        region.expand(amount);
        return *this;
    }

    /**
     * @brief Returns a copy of this GenomicRegion expanded by the specified amount.
     *
     * Creates a new GenomicRegion that is expanded by the provided amount
     * without modifying the original GenomicRegion.
     *
     * @param amount The amount to expand the region
     * @return A new GenomicRegion object with expanded boundaries
     */
    [[nodiscard]] auto expanded(const size_t amount) const -> GenomicRegion {
        auto expandedRegion = *this;
        expandedRegion.expand(amount);
        return expandedRegion;
    }

    // Operators
    auto operator<(const GenomicRegion& other) const -> bool {
        if (referenceIDIndex != other.referenceIDIndex) {
            return referenceIDIndex < other.referenceIDIndex;
        }
        return region.startPosition < other.region.startPosition;
    }

    auto operator>(const GenomicRegion& other) const -> bool { return other < *this; }

    auto operator==(const GenomicRegion& other) const -> bool {
        return referenceIDIndex == other.referenceIDIndex && region == other.region &&
               strand == other.strand;
    }

    /**
     * @brief Checks whether this GenomicRegion overlaps with another GenomicRegion, considering
     * orientation and tolerance.
     *
     * This function first verifies if both GenomicRegions refer to the same reference ID. If they
     * do, it checks the specified orientation and uses the region.overlaps method with the provided
     * tolerance to determine if there is any overlap.
     *
     * @param other The other GenomicRegion to compare with.
     * @param orientation The orientation to consider when checking for overlaps (SAME, OPPOSITE, or
     * BOTH).
     * @param tolerance The tolerance value to apply when checking for overlaps in the underlying
     * region.
     * @return True if the GenomicRegions overlap under the specified orientation and tolerance;
     * false otherwise.
     */
    [[nodiscard]] auto overlapsWithTolerance(const GenomicRegion& other,
                                             const dataTypes::GenomicOrientation orientation,
                                             int tolerance) const noexcept -> bool {
        if (referenceIDIndex != other.referenceIDIndex) {
            return false;
        }

        switch (orientation) {
            case dataTypes::GenomicOrientation::SAME:
                return strand == other.strand &&
                       region.overlapsWithTolerance(other.region, tolerance);
            case dataTypes::GenomicOrientation::OPPOSITE:
                return strand == !other.strand &&
                       region.overlapsWithTolerance(other.region, tolerance);
            case dataTypes::GenomicOrientation::BOTH:
                return region.overlapsWithTolerance(other.region, tolerance);
            default:
                return false;
        }
    };

    [[nodiscard]] auto overlapsWithShortestSegmentFraction(
        const GenomicRegion& other, const dataTypes::GenomicOrientation orientation,
        const float shortestOverlapFraction) const noexcept -> bool {
        if (referenceIDIndex != other.referenceIDIndex) {
            return false;
        }

        switch (orientation) {
            case dataTypes::GenomicOrientation::SAME:
                return strand == other.strand && region.overlapsWithShortestRegionFraction(
                                                     other.region, shortestOverlapFraction);
            case dataTypes::GenomicOrientation::OPPOSITE:
                return strand == !other.strand && region.overlapsWithShortestRegionFraction(
                                                      other.region, shortestOverlapFraction);
            case dataTypes::GenomicOrientation::BOTH:
                return region.overlapsWithShortestRegionFraction(other.region,
                                                                 shortestOverlapFraction);
            default:
                return false;
        }
    }

    /**
     * @brief Merges this GenomicRegion with another GenomicRegion if they overlap.
     *
     * This method attempts to merge the current GenomicRegion with the specified 'other'
     * GenomicRegion. The merge is performed only if the two regions overlap under the
     * annotation::Orientation::SAME orientation, while considering a tolerance of 1 for
     * blunt-end boundaries. If the regions do not overlap, no merge is performed and false
     * is returned.
     *
     * @param other The GenomicRegion to attempt merging with.
     * @return True if the regions were merged successfully, false otherwise.
     */
    [[nodiscard]] auto merge(const GenomicRegion& other,
                             const GenomicStrandSpecificity& strandSpecificity =
                                 GenomicStrandSpecificity::SPECIFIC) noexcept -> bool {
        // Tolerance 1 because blunt ended regions can be merged
        if (!overlapsWithTolerance(
                other, GenomicOrientation::fromStrandSpecificity(strandSpecificity), 1)) {
            return false;
        }

        if (strand != other.strand) {
            strand = GenomicStrand::NONE;
        }

        region.merge(other.region);
        return true;
    }

    /**
     * @brief Combines this GenomicRegion with another if they share the same strand and reference
     * ID.
     *
     * If the strands and reference IDs match, the regions are merged and the method returns true.
     * Otherwise, no changes are made and the method returns false.
     *
     * @param other The GenomicRegion to combine with.
     * @return True if the region was updated; false otherwise.
     */
    [[nodiscard]] auto combine(const GenomicRegion& other) noexcept -> bool {
        if (referenceIDIndex != other.referenceIDIndex ||
            (strand != GenomicStrand::NONE && strand != other.strand)) {
            return false;
        }

        const auto prevRegion = region;

        region.merge(other.region);
        return prevRegion != region;
    }

    /**
     * @brief Computes the fraction of overlap between this GenomicRegion and another.
     *
     * This function first checks if the two regions overlap under the specified orientation.
     * If they do not overlap, 0.0 is returned. Otherwise, the function calculates the overlapping
     * length and divides it by the length of the current region. This approach avoids repeated
     * calls to getters to improve efficiency.
     *
     * @param other The other GenomicRegion to compare with.
     * @param orientation The specified orientation for overlap checking.
     * @return The fraction of overlap for this region as a double.
     */
    [[nodiscard]] auto overlapFraction(
        const GenomicRegion& other, const dataTypes::GenomicOrientation orientation) const noexcept
        -> float {
        // Overlaps at least one base
        if (!overlapsWithTolerance(other, orientation, 0)) {
            return 0.0;
        }

        const auto start = getStart();
        const auto end = getEnd();
        const auto otherStart = other.getStart();
        const auto otherEnd = other.getEnd();

        const auto overlapStart = std::max(start, otherStart);
        const auto overlapEnd = std::min(end, otherEnd);
        const auto overlapLength = overlapEnd - overlapStart;

        return static_cast<float>(overlapLength) / static_cast<float>(length());
    }

    [[nodiscard]] auto contains(const int32_t position) const noexcept -> bool {
        return region.contains(position);
    }

    [[nodiscard]] auto contains(const GenomicRegion& other,
                                const dataTypes::GenomicOrientation orientation) const noexcept
        -> bool {
        if (referenceIDIndex != other.referenceIDIndex) {
            return false;
        }

        switch (orientation) {
            case dataTypes::GenomicOrientation::SAME:
                return strand == other.strand && region.contains(other.region.startPosition) &&
                       region.contains(other.region.endPosition);
            case dataTypes::GenomicOrientation::OPPOSITE:
                return strand == !other.strand && region.contains(other.region.startPosition) &&
                       region.contains(other.region.endPosition);
            case dataTypes::GenomicOrientation::BOTH:
                return region.contains(other.region.startPosition) &&
                       region.contains(other.region.endPosition);
            default:
                return false;
        }
    }

   private:
    int referenceIDIndex;
    Region region;
    GenomicStrand strand;
};

inline auto operator<<(std::ostream& outputStream, const GenomicRegion& genomicRegion)
    -> std::ostream& {
    return outputStream << genomicRegion.getReferenceIDIndex() << ":" << genomicRegion.getStart()
                        << "-" << genomicRegion.getEnd() << ' ' << genomicRegion.getStrand()
                        << '\n';
};
}  // namespace dataTypes
