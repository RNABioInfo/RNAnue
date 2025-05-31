#pragma once

// Standard
#include <cstddef>
#include <vector>

// Internal
#include "ReadAlignment.hpp"
#include "SplitRecords.hpp"

namespace dataTypes {

struct ReadAlignments {
    std::vector<ReadAlignmentVariant> alignments;

    [[nodiscard]] constexpr auto alignmentCount() const noexcept -> size_t {
        return alignments.size();
    };
};

using SplitRecordsVariantGroups = std::vector<MultimericRecords>;
}  // namespace dataTypes
