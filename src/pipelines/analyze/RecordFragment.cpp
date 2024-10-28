#include "RecordFragment.hpp"

// Standard
#include <optional>

// Include
#include "CustomSamTags.hpp"
#include "GenomicRegion.hpp"
#include "Utility.hpp"

namespace pipelines::analyze {

auto RecordFragment::fromSamRecord(const SamRecord &record) -> std::optional<RecordFragment> {
    if (!record.reference_position().has_value() || !record.reference_id().has_value()) {
        return std::nullopt;
    }

    const auto isReverseStrand =
        static_cast<bool>(record.flag() & seqan3::sam_flag::on_reverse_strand);
    const dataTypes::Strand strand{isReverseStrand ? dataTypes::Strand::REVERSE
                                                   : dataTypes::Strand::FORWARD};

    const auto start = record.reference_position();
    const std::optional<int32_t> end = dataTypes::recordEndPosition(record);

    if (!start.has_value() || !end.has_value()) {
        return std::nullopt;
    }

    const double hybridizationEnergy = record.tags().get<"XE"_tag>();
    const double complementarityScore = record.tags().get<"XC"_tag>();

    return RecordFragment{.recordID = record.id(),
                          .referenceIDIndex = record.reference_id().value(),
                          .strand = strand,
                          .start = start.value(),
                          .end = end.value(),
                          .complementarityScore = complementarityScore,
                          .hybridizationEnergy = hybridizationEnergy};
}

auto RecordFragment::toGenomicRegion(const std::deque<std::string> &referenceIDs) const
    -> dataTypes::GenomicRegion {
    return dataTypes::GenomicRegion{referenceIDs[referenceIDIndex], start, end, strand};
}

auto RecordFragment::toFeature(const std::deque<std::string> &referenceIDs,
                               const std::string &featureID, const std::string &featureType) const
    -> dataTypes::GenomicFeature {
    return dataTypes::GenomicFeature{.referenceID = referenceIDs[referenceIDIndex],
                                     .type = featureType,
                                     .startPosition = start,
                                     .endPosition = end,
                                     .strand = strand,
                                     .id = featureID,
                                     .groupID = std::nullopt,
                                     .geneName = std::nullopt};
}

void RecordFragment::merge(const RecordFragment &other) {
    start = std::min(start, other.start);
    end = std::max(end, other.end);
    complementarityScore = std::max(complementarityScore, other.complementarityScore);
    hybridizationEnergy = std::min(hybridizationEnergy, other.hybridizationEnergy);
}

auto RecordFragment::operator==(const RecordFragment &other) const -> bool {
    constexpr double EPSILON = 1e-6;
    return referenceIDIndex == other.referenceIDIndex && strand == other.strand &&
           start == other.start && end == other.end &&
           helper::isEqual(complementarityScore, other.complementarityScore, EPSILON) &&
           helper::isEqual(hybridizationEnergy, other.hybridizationEnergy, EPSILON);
}

auto RecordFragment::operator<(const RecordFragment &other) const -> bool {
    if (referenceIDIndex != other.referenceIDIndex) {
        return referenceIDIndex < other.referenceIDIndex;
    }

    return start < other.start;
}

auto RecordFragment::operator>(const RecordFragment &other) const -> bool { return other < *this; }

}  // namespace pipelines::analyze
