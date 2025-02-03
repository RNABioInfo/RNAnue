#pragma once

// Standard
#include <cstdint>
#include <memory>
#include <variant>

// Internal
#include "FeatureAnnotator.hpp"
#include "GenomicOrientation.hpp"

namespace SplitRecordsEvaluationParameters {
struct BaseParameters {
    double minComplementarity;
    double minComplementarityFraction;
    double mfeThreshold;
    bool includeWobbleBasePairsInCrosslinkingSites;
};

struct SplicingParameters {
    BaseParameters baseParameters;
    dataTypes::GenomicOrientation orientation;
    int32_t splicingTolerance;
    bool allowAlternativeSplicing;
    std::shared_ptr<const annotation::FeatureAnnotator> featureAnnotator;
};

using ParameterVariant = std::variant<SplitRecordsEvaluationParameters::BaseParameters,
                                      SplitRecordsEvaluationParameters::SplicingParameters>;
}  // namespace SplitRecordsEvaluationParameters
