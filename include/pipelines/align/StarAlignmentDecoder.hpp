#pragma once

#include <array>
#include <deque>
#include <vector>
#include "SamRecord.hpp"

namespace pipelines::align {
struct StarAlignmentHypothesis {
    int hitIndex{};
    bool paired{};
    std::vector<size_t> rawIndices;
    std::vector<dataTypes::SamRecord> fragments;
    std::array<size_t, 2> queryLengths{};
    std::array<size_t, 2> coveredBases{};
};

// Shared by output filtering and detection. No mapper tags are fabricated.
class StarAlignmentDecoder {
   public:
    explicit StarAlignmentDecoder(std::deque<std::string> referenceIDs,
                                  std::vector<size_t> referenceLengths = {})
        : referenceIDs(std::move(referenceIDs)), referenceLengths(std::move(referenceLengths)) {}
    [[nodiscard]] auto decode(const std::vector<dataTypes::SamRecord>& records,
                              bool excludeSoftClipping) const
        -> std::vector<StarAlignmentHypothesis>;

   private:
    std::deque<std::string> referenceIDs;
    std::vector<size_t> referenceLengths;
};
}  // namespace pipelines::align
