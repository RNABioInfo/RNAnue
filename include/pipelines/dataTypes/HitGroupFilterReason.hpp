#pragma once

// Standard
#include <cstdint>
#include <ostream>

namespace dataTypes {

/**
 * @enum HitGroupFilterReason
 * @brief Enumerates the various reasons for a split hit group filtering.
 *
 * Each enumerator represents a specific type of failure case encountered in processing hit groups.
 * Other then HitGroupFailureReason this results in the downgrade of the hit group from split to
 * multiple singleton HitGroups.
 */
enum class HitGroupFilterReason : std::uint8_t {
    SPLICING,
    COMPLEMENTARITY,
    HYBRIDIZATION,
    CROSSLINKING_COUNT
};

inline auto operator<<(std::ostream& ostream, const HitGroupFilterReason& reason) -> std::ostream& {
    switch (reason) {
        case HitGroupFilterReason::SPLICING:
            ostream << "Splicing";
            break;
        case HitGroupFilterReason::COMPLEMENTARITY:
            ostream << "Complementarity";
            break;
        case HitGroupFilterReason::HYBRIDIZATION:
            ostream << "Hybridization";
            break;
        case HitGroupFilterReason::CROSSLINKING_COUNT:
            ostream << "Crosslinking Count";
            break;
    }
    return ostream;
};

}  // namespace dataTypes
