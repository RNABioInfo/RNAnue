#include "GenomicRegion.hpp"

// Standard
#include <optional>

// seqan3
#include <seqan3/io/sam_file/sam_flag.hpp>

// Internal
#include "GenomicStrand.hpp"
#include "SamRecord.hpp"

namespace dataTypes {

auto GenomicRegion::fromSamRecord(const SamRecord &record) -> std::optional<GenomicRegion> {
    const auto start = record.reference_position();
    const auto end = recordEndPosition(record);

    if (!start.has_value() || !end.has_value() || !record.reference_id()) {
        return std::nullopt;
    }

    const auto isReverseStrand =
        static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand);
    const GenomicStrand strand{isReverseStrand ? GenomicStrand::REVERSE : GenomicStrand::FORWARD};

    return GenomicRegion{
        record.reference_id().value(), {.startPosition = *start, .endPosition = *end}, strand};
}

}  // namespace dataTypes
