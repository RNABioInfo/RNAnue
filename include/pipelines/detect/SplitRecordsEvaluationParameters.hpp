#pragma once

// Standard
#include <cstdint>
#include <variant>

// Internal
#include "FeatureAnnotator.hpp"
#include "Orientation.hpp"

namespace SplitRecordsEvaluationParameters {
struct BaseParameters {
    double minComplementarity;
    double minComplementarityFraction;
    double mfeThreshold;
    bool includeWobbleBasePairsInCrosslinkingSites;
};

struct SplicingParameters {
    BaseParameters baseParameters;
    annotation::Orientation orientation;
    int32_t splicingTolerance;
    const annotation::FeatureAnnotator* featureAnnotator;
};

using ParameterVariant = std::variant<SplitRecordsEvaluationParameters::BaseParameters,
                                      SplitRecordsEvaluationParameters::SplicingParameters>;
}  // namespace SplitRecordsEvaluationParameters
