#pragma once

// Standard
#include <cstdint>
#include <optional>
#include <ostream>
#include <string>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "SamRecord.hpp"

namespace dataTypes {

using seqan3::operator""_tag;

struct RecordFragment {
    std::string recordID;
    GenomicRegion genomicRegion;
    double complementarityScore;
    double hybridizationEnergy;
    int32_t interCrosslinkingSiteCount;
    float transcriptContribution{1.0F};

    [[nodiscard]] static auto fromSamRecord(const SamRecord &record)
        -> std::optional<RecordFragment>;

    [[nodiscard]] auto toFeature(const std::string &featureID,
                                 const std::string &featureType = "transcript") const
        -> dataTypes::GenomicFeature;

    [[nodiscard]] auto merge(const RecordFragment &other) -> bool;

    auto operator==(const RecordFragment &other) const -> bool;

    auto operator<(const RecordFragment &other) const -> bool;
    auto operator>(const RecordFragment &other) const -> bool;
};

inline auto operator<<(std::ostream &ostream, const RecordFragment &segment) -> std::ostream & {
    return ostream << "Segment{recordID: " << segment.recordID
                   << "\n, genomicRegion: " << segment.genomicRegion
                   << ", maxComplementarityScore: " << segment.complementarityScore
                   << "\n, minHybridizationEnergy: " << segment.hybridizationEnergy << "}\n";
};

}  // namespace dataTypes
