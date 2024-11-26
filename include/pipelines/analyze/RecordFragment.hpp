#pragma once

// Standard
#include <cstdint>
#include <deque>
#include <optional>
#include <string>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "SamRecord.hpp"

using seqan3::operator""_tag;
using namespace dataTypes;

namespace pipelines::analyze {

struct RecordFragment {
    std::string recordID;
    int32_t referenceIDIndex;
    dataTypes::GenomicStrand strand;
    int32_t start;
    int32_t end;
    double complementarityScore;
    double hybridizationEnergy;
    int32_t crosslinkingSiteCount;

    static auto fromSamRecord(const SamRecord &record) -> std::optional<RecordFragment>;

    [[nodiscard]] auto toGenomicRegion(const std::deque<std::string> &referenceIDs) const
        -> dataTypes::GenomicRegion;

    [[nodiscard]] auto toFeature(const std::deque<std::string> &referenceIDs,
                                 const std::string &featureID,
                                 const std::string &featureType = "transcript") const
        -> dataTypes::GenomicFeature;

    void merge(const RecordFragment &other);

    auto operator==(const RecordFragment &other) const -> bool;

    auto operator<(const RecordFragment &other) const -> bool;
    auto operator>(const RecordFragment &other) const -> bool;
};

inline auto operator<<(std::ostream &ostream, const RecordFragment &segment) -> std::ostream & {
    return ostream << "Segment{recordID: " << segment.recordID
                   << ", referenceIDIndex: " << segment.referenceIDIndex
                   << ", strand: " << segment.strand << ", start: " << segment.start
                   << ", end: " << segment.end
                   << ", maxComplementarityScore: " << segment.complementarityScore
                   << ", minHybridizationEnergy: " << segment.hybridizationEnergy << "}";
};

}  // namespace pipelines::analyze
