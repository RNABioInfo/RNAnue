#pragma once

#include <algorithm>

#include "GenomicRegion.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "Logger.hpp"

namespace dataTypes {

struct SortedGenomicRegionPair {
    /**
     * @brief Constructs a SortedGenomicRegionPair by sorting the input regions.
     *
     * This constructor ensures that the first region starts before the second region.
     *
     * @param regionOne The first region to be sorted
     * @param regionTwo The second region to be sorted
     */
    constexpr SortedGenomicRegionPair(GenomicRegion regionOne, GenomicRegion regionTwo)
        : firstRegion(std::min(regionOne, regionTwo)),
          secondRegion(std::max(regionOne, regionTwo)) {}

    GenomicRegion firstRegion;
    GenomicRegion secondRegion;

    auto operator==(const SortedGenomicRegionPair& other) const noexcept -> bool {
        return firstRegion == other.firstRegion && secondRegion == other.secondRegion;
    }

    /**
     * @brief Merges this SortedGenomicRegionPair with another one.
     *
     * Tries to merge the first and second GenomicRegion members of this object
     * with the corresponding members of the given pair. If either merge operation
     * fails, the changes are rolled back and the method returns false.
     *
     * @param other The other SortedGenomicRegionPair to be merged.
     * @return true if both merges succeed, false otherwise.
     */
    auto merge(const SortedGenomicRegionPair& other,
               const GenomicStrandSpecificity& strandSpecificity) noexcept -> bool {
        const auto tmpFirstRegion = firstRegion;
        const auto tmpSecondRegion = secondRegion;

        if ((!firstRegion.merge(other.firstRegion, strandSpecificity) ||
             !secondRegion.merge(other.secondRegion, strandSpecificity))) [[unlikely]] {
            Logger::log<LogLevel::WARNING>(
                "Attempted to merge non-overlapping region pairs with mode: ", strandSpecificity,
                "\n", firstRegion, secondRegion, other.firstRegion, other.secondRegion);

            firstRegion = tmpFirstRegion;
            secondRegion = tmpSecondRegion;

            return false;
        }

        return true;
    }
};

}  // namespace dataTypes
