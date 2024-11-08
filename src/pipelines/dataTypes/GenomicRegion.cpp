#include "GenomicRegion.hpp"

namespace dataTypes {

auto GenomicRegion::fromSamRecord(const dataTypes::SamRecord &record,
                                  const std::deque<std::string> &referenceIDs)
    -> std::optional<GenomicRegion> {
    const auto start = record.reference_position();
    const auto end = recordEndPosition(record);

    if (!start.has_value() || !end.has_value()) {
        return std::nullopt;
    }

    const auto isReverseStrand =
        static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand);
    const GenomicStrand strand{isReverseStrand ? GenomicStrand::REVERSE : GenomicStrand::FORWARD};

    return GenomicRegion{referenceIDs[record.reference_id().value()], start.value(),
                         end.value() + 1, strand};
}

auto GenomicRegion::fromGenomicFeature(const dataTypes::GenomicFeature &feature) -> GenomicRegion {
    return GenomicRegion{feature.referenceID, feature.startPosition, feature.endPosition,
                         feature.strand};
}

}  // namespace dataTypes
