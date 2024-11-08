#pragma once

// Standard
#include <deque>
#include <optional>
#include <string>
#include <utility>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicStrand.hpp"
#include "SamRecord.hpp"

namespace dataTypes {
/**
 * @brief Represents a genomic region.
 */
struct GenomicRegion {
    std::string referenceID;
    int32_t startPosition;
    int32_t endPosition;
    std::optional<GenomicStrand> strand;

    /**
     * @brief Constructs a GenomicRegion object.
     * @param referenceID The reference ID of the genomic region.
     * @param startPosition The start position of the genomic region.
     * @param endPosition The end position of the genomic region (exclusive).
     * @param strand The strand of the genomic region (optional).
     */
    GenomicRegion(std::string referenceID, int32_t startPosition, int32_t endPosition,
                  std::optional<GenomicStrand> strand = std::nullopt)
        : referenceID(std::move(referenceID)),
          startPosition(startPosition),
          endPosition(endPosition),
          strand(strand) {};

    /**
     * @brief Creates a GenomicRegion object from a SamRecord.
     * @param record The SamRecord object.
     * @param referenceIDs The deque of reference IDs.
     * @return An optional GenomicRegion object.
     */
    static auto fromSamRecord(const dataTypes::SamRecord &record,
                              const std::deque<std::string> &referenceIDs)
        -> std::optional<GenomicRegion>;

    static auto fromGenomicFeature(const dataTypes::GenomicFeature &feature) -> GenomicRegion;
};

inline auto operator<<(std::ostream &outputStream, const GenomicRegion &region) -> std::ostream & {
    return outputStream << region.referenceID << ":" << region.startPosition << "-"
                        << region.endPosition << ' '
                        << (region.strand.has_value() ? region.strand.value()
                                                      : GenomicStrand::FORWARD)
                        << '\n';
};
}  // namespace dataTypes
