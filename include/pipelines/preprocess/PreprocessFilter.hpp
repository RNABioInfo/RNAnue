#pragma once

// standard
#include <cstddef>

// Internal
#include "FastqRecord.hpp"
#include "SequenceQualityAlgorithms.hpp"

namespace pipelines::preprocess::PreprocessFilter {

using namespace dataTypes;

struct Criteria {
    size_t minLengthThreshold;
    size_t minQualityThreshold;
};

static inline auto passes(const Criteria &criteria, const FastqRecord &record) -> bool {
    const auto meanQual = SequenceQualityAlgorithms::meanQualityScore(record.base_qualities());

    const bool passesQual = meanQual >= double(criteria.minQualityThreshold);

    const bool passesLen = record.sequence().size() >= criteria.minLengthThreshold;

    return passesQual && passesLen;
}

}  // namespace pipelines::preprocess::PreprocessFilter
