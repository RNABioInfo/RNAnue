#pragma once

// Standard
#include <algorithm>
#include <vector>

// Internal
#include "SamRecord.hpp"

namespace dataTypes {

// Represents a single chimeric read which is split up into its blocks
struct SplitRecords : public std::vector<SamRecord> {
    SplitRecords(const std::vector<SamRecord>& records) : std::vector<SamRecord>(records) {
        std::sort(begin(), end(), dataTypes::operator<);
    }

    SplitRecords() = default;

    auto operator<(const SplitRecords& rhs) const -> bool {
        const auto& lhs_ref_id = this->back().reference_id();
        const auto& rhs_ref_id = rhs.back().reference_id();

        if (lhs_ref_id != rhs_ref_id) {
            return lhs_ref_id < rhs_ref_id;
        }

        const auto lhsEnd = recordEndPosition(this->back());
        const auto rhsEnd = recordEndPosition(rhs.back());

        return lhsEnd < rhsEnd;
    }

    auto operator>(const SplitRecords& rhs) const -> bool { return rhs < *this; }
};

// Split records are sorted by the end position of the last split record

}  // namespace dataTypes
