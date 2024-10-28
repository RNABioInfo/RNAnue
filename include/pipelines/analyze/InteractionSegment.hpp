#pragma once

// Standard
#include <cstdint>
#include <ostream>

// Internal
#include "GenomicStrand.hpp"
#include "RecordFragment.hpp"

namespace pipelines::analyze {

struct InteractionSegment {
    InteractionSegment(const RecordFragment& recordFragment)
        : referenceIDIndex(recordFragment.referenceIDIndex),
          strand(recordFragment.strand),
          start(recordFragment.start),
          end(recordFragment.end) {}

    InteractionSegment(int32_t referenceIDIndex, dataTypes::Strand strand, int32_t start,
                       int32_t end)
        : referenceIDIndex(referenceIDIndex), strand(strand), start(start), end(end) {}

    // Getters
    [[nodiscard]] auto getReferenceIDIndex() const { return referenceIDIndex; }
    [[nodiscard]] auto getStrand() const { return strand; }
    [[nodiscard]] auto getStart() const { return start; }
    [[nodiscard]] auto getEnd() const { return end; }

    auto operator<(const InteractionSegment& other) const -> bool {
        if (referenceIDIndex != other.referenceIDIndex) {
            return referenceIDIndex < other.referenceIDIndex;
        }
        return start < other.start;
    }

    auto operator>(const InteractionSegment& other) const -> bool { return other < *this; }

    auto operator==(const InteractionSegment& other) const -> bool {
        return referenceIDIndex == other.referenceIDIndex && start == other.start &&
               end == other.end && strand == other.strand;
    }

    [[nodiscard]] auto overlaps(const InteractionSegment& other, int graceDistance) const -> bool {
        if (referenceIDIndex != other.referenceIDIndex || strand != other.strand) {
            return false;
        }

        return (start <= other.start) ? (other.start - end <= graceDistance)
                                      : (start - other.end <= graceDistance);
    };

    void merge(const InteractionSegment& other) {
        start = std::min(start, other.start);
        end = std::max(end, other.end);
    }

   private:
    int32_t referenceIDIndex;
    dataTypes::Strand strand;
    int32_t start;
    int32_t end;
};

inline auto operator<<(std::ostream& outputStream, const InteractionSegment& segment)
    -> std::ostream& {
    outputStream << "InteractionSegment(referenceIDIndex: " << segment.getReferenceIDIndex()
                 << ", strand: " << segment.getStrand() << ", start: " << segment.getStart()
                 << ", end: " << segment.getEnd() << ")";
    return outputStream;
}

}  // namespace pipelines::analyze
