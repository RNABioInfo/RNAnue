#pragma once

// Standard
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <variant>

// Internal
#include "DetectParameters.hpp"
#include "FeatureAnnotator.hpp"
#include "GenomicOrientation.hpp"
#include "OneOf.hpp"

namespace pipelines::detect::ReadGroupEvaluationParameters {

struct Base {
    size_t maxPrimariyAlignmentsCount;
    double minComplementarity;
    double minComplementarityFraction;
    double mfeThreshold;
    bool includeWobbleBasePairsInCrosslinkingSites;
    double minHitGroupContribution;
    std::shared_ptr<const annotation::FeatureAnnotator> featureAnnotator;
    annotation::GenomicOrientation featureOrientation;
};

struct Splicing {
    Splicing(Base baseParams, std::int32_t splicingTolerance,
             bool allowAlternativeSplicing) noexcept
        : baseParameters(std::move(baseParams)),
          splicingTolerance(splicingTolerance),
          allowAlternativeSplicing(allowAlternativeSplicing) {}

    Base baseParameters;
    int32_t splicingTolerance;
    bool allowAlternativeSplicing;
};

using ParameterVariant = std::variant<Base, Splicing>;

[[nodiscard]] inline auto makeParams(
    const DetectParameters& params,
    std::shared_ptr<const annotation::FeatureAnnotator> featureAnnotator) -> ParameterVariant {
    ReadGroupEvaluationParameters::Base baseParams = {
        .maxPrimariyAlignmentsCount = params.maxPrimaryAlignmentCount,
        .minComplementarity = params.minimumComplementarity,
        .minComplementarityFraction = params.minimumSiteLengthRatio,
        .mfeThreshold = params.maxHybridizationEnergy,
        .includeWobbleBasePairsInCrosslinkingSites =
            params.includeWobbleBasePairsInCrosslinkingSites,
        .minHitGroupContribution = params.minHitGroupContribution,
        .featureAnnotator = std::move(featureAnnotator),
        .featureOrientation = params.featureOrientation};

    if (!params.removeSplicingEvents) {
        return baseParams;
    }

    return ReadGroupEvaluationParameters::Splicing{baseParams, params.splicingTolerance,
                                                   params.removeAlternativeSplicing};
}

template <typename T>
concept Type = one_of<T, Base, Splicing>;

}  // namespace pipelines::detect::ReadGroupEvaluationParameters
