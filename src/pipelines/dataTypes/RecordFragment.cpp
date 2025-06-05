#include "RecordFragment.hpp"

// Standard
#include <algorithm>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>

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
                          .transcriptContribution = transcriptContribution};
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
