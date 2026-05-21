#include "RecordFragment.hpp"

// Standard
#include <algorithm>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

// Include
#include "CustomSamTags.hpp"  // IWYU pragma: keep
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "SamRecord.hpp"
#include "Utility.hpp"

namespace dataTypes {
namespace {

[[nodiscard]] auto alignedCoverageIntervals(const SamRecord &record)
    -> std::vector<pipelines::analyze::CoverageInterval> {
    std::vector<pipelines::analyze::CoverageInterval> intervals;

    const auto start = record.reference_position();
    if (!start.has_value()) {
        return intervals;
    }

    int32_t referencePosition = *start;
    for (const auto &cigar : record.cigar_sequence()) {
        const auto length = static_cast<int32_t>(get<0>(cigar));

        if (cigar == 'M'_cigar_operation || cigar == '='_cigar_operation ||
            cigar == 'X'_cigar_operation) {
            intervals.push_back(
                {.start = referencePosition, .end = referencePosition + length});
            referencePosition += length;
        } else if (cigar == 'D'_cigar_operation || cigar == 'N'_cigar_operation) {
            referencePosition += length;
        }
    }

    return intervals;
}

}  // namespace

auto RecordFragment::fromSamRecord(const SamRecord &record) -> std::optional<RecordFragment> {
    const std::optional<GenomicRegion> genomicRegion = GenomicRegion::fromSamRecord(record);

    if (!genomicRegion) {
        return std::nullopt;
    }

    const double hybridizationEnergy = record.tags().get<"XE"_tag>();
    const double complementarityScore = record.tags().get<"XC"_tag>();
    const int32_t interCrosslinkingSiteCount = record.tags().get<"XO"_tag>();
    const float transcriptContribution = record.tags().get<"XB"_tag>();

    return RecordFragment{.recordID = record.id(),
                          .genomicRegion = *genomicRegion,
                          .complementarityScore = complementarityScore,
                          .hybridizationEnergy = hybridizationEnergy,
                          .interCrosslinkingSiteCount = interCrosslinkingSiteCount,
                          .transcriptContribution = transcriptContribution,
                          .coverageIntervals = alignedCoverageIntervals(record)};
}

auto RecordFragment::toFeature(const std::string &featureID, const std::string &featureType) const
    -> GenomicFeature {
    return GenomicFeature{featureType, genomicRegion, featureID, std::nullopt, std::nullopt};
}

[[nodiscard]] auto RecordFragment::merge(const RecordFragment &other) -> bool {
    if (!genomicRegion.merge(other.genomicRegion)) {
        return false;
    }

    complementarityScore = std::max(complementarityScore, other.complementarityScore);
    hybridizationEnergy = std::min(hybridizationEnergy, other.hybridizationEnergy);
    coverageIntervals.insert(coverageIntervals.end(), other.coverageIntervals.begin(),
                             other.coverageIntervals.end());

    return true;
}

auto RecordFragment::operator==(const RecordFragment &other) const -> bool {
    constexpr double EPSILON = 1e-6;
    return genomicRegion == other.genomicRegion &&
           helper::isApproxEqual(complementarityScore, other.complementarityScore, EPSILON) &&
           helper::isApproxEqual(hybridizationEnergy, other.hybridizationEnergy, EPSILON);
}

auto RecordFragment::operator<(const RecordFragment &other) const -> bool {
    return genomicRegion < other.genomicRegion;
}

auto RecordFragment::operator>(const RecordFragment &other) const -> bool { return other < *this; }

}  // namespace dataTypes
