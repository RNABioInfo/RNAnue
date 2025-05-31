#pragma once

// Standard
#include <cstdint>
#include <ostream>

namespace dataTypes {

/**
 * @enum HitGroupFailureReason
 * @brief Enumerates the various reasons for a hit group failure.
 *
 * Each enumerator represents a specific type of failure case encountered in processing hit groups.
 * Other then HitGroupFilterReason this results in the removal of the hit group from consideration.
 */
enum class HitGroupFailureReason : std::uint8_t {
    MALFORMED_RECORD,
    MALFORMED_HIT_GROUP,
    UNMAPPED,
    MULTIMAPPING,
    MAPPING_QUALITY,
    FRAGMENT_LENGTH,
    MIN_CONTRIBUTION,
    NOT_IMPLEMENTED,
    FAILED_COMPLEMENTARITY,
    FAILED_HYBRIDIZATION
};

inline auto operator<<(std::ostream& ostream, const HitGroupFailureReason& reason)
    -> std::ostream& {
    switch (reason) {
        case HitGroupFailureReason::MALFORMED_RECORD:
            ostream << "Malformed Record";
            break;
        case HitGroupFailureReason::MALFORMED_HIT_GROUP:
            ostream << "Malformed Hit Group";
            break;
        case HitGroupFailureReason::UNMAPPED:
            ostream << "Unmapped";
            break;
        case HitGroupFailureReason::MULTIMAPPING:
            ostream << "Multimapping";
            break;
        case HitGroupFailureReason::MAPPING_QUALITY:
            ostream << "Mapping Quality";
            break;
        case HitGroupFailureReason::FRAGMENT_LENGTH:
            ostream << "Fragment Length";
            break;
        case HitGroupFailureReason::MIN_CONTRIBUTION:
            ostream << "Minimum Contribution";
            break;
        case HitGroupFailureReason::NOT_IMPLEMENTED:
            ostream << "Not Implemented";
            break;
        case HitGroupFailureReason::FAILED_COMPLEMENTARITY:
            ostream << "Failed Complementarity";
            break;
        case HitGroupFailureReason::FAILED_HYBRIDIZATION:
            ostream << "Failed Hybridization";
            break;
    }
    return ostream;
};

}  // namespace dataTypes
