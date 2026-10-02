#pragma once

// Standard
#include <cfloat>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>
#include <tuple>

// Internal
#include "ArithmeticParameterOption.hpp"
#include "DefaultedParameterOption.hpp"

struct DetectOptions {
    static constexpr std::string_view optionsDescription = "Detect Pipeline"sv;

    static constexpr ArithmeticParameterOption<size_t, 5, {.lowerBound = 1, .upperBound = SIZE_MAX}>
        maxPrimaryAlignmentCount{{.shortName = std::nullopt, .longName = "max_primary_alignments"},
                                 "maximum count of primary alignments per read"sv};

    static constexpr ArithmeticParameterOption<size_t, 0, {.lowerBound = 0, .upperBound = 100}>
        minMappingQuality{{.shortName = std::nullopt, .longName = "min_mapping_quality"},
                          "minimum mapping quality for reads to be considered in the analysis"sv};

    static constexpr ArithmeticParameterOption<double, 0.5, {.lowerBound = 0.0, .upperBound = 1.0}>
        minComplementarity{{.shortName = std::nullopt, .longName = "min_complementarity"},
                           "complementarity cutoff for split reads"sv};

    static constexpr ArithmeticParameterOption<double, 0.1, {.lowerBound = 0.0, .upperBound = 1.0}>
        siteLengthRatio{{.shortName = std::nullopt, .longName = "min_site_length_ratio"},
                        "fraction of shortest record match cutoff"sv};

    static constexpr ArithmeticParameterOption<size_t, 20,
                                               {.lowerBound = 0, .upperBound = SIZE_MAX}>
        minDetectLength{{.shortName = std::nullopt, .longName = "min_detect_length"},
                        "minimum fragment length after clipping"sv};

    static constexpr DefaultedParameterOption<double, 0.0> maxEnergy{
        {.shortName = std::nullopt, .longName = "max_hybridization_energy"},
        "hybridization energy cutoff for split reads"sv};

    static constexpr DefaultedParameterOption<bool, false> excludeSoftClipping{
        {.shortName = std::nullopt, .longName = "exclude_soft_clipping"},
        "exclude soft clipping from the alignments"sv};

    static constexpr DefaultedParameterOption<bool, false> filterSplicing{
        {.shortName = std::nullopt, .longName = "filter_splicing"},
        "splicing events are removed in the detection of split reads"sv};

    static constexpr DefaultedParameterOption<bool, true> removeAltSplicing{
        {.shortName = std::nullopt,
         .longName = "remove_alt_splicing"},
        "remove alternative splicing events"sv};

    static constexpr DefaultedParameterOption<int, 5> splicingTolerance{
        {.shortName = std::nullopt, .longName = "splicing_tolerance"}, "tolerance for splicing events"sv};

    static constexpr DefaultedParameterOption<bool, false> includeWobble{
        {.shortName = std::nullopt, .longName = "include_wobble"},
        "[EXPERIMENTAL] wobble base pairs are allowed in crosslinking site evaluation"sv};

    static constexpr ArithmeticParameterOption<double, 0.1, {.lowerBound = 0, .upperBound = 1}>
        minHitGroupContribution{
            {.shortName = std::nullopt, .longName = "min_hit_contribution"},
            "minimum contribution of alignment within all alignments of this read"sv};

    static constexpr auto allOptions = std::make_tuple(
        maxPrimaryAlignmentCount, minMappingQuality, minComplementarity, siteLengthRatio,
        minDetectLength, maxEnergy, excludeSoftClipping, filterSplicing, removeAltSplicing,
        splicingTolerance, includeWobble, minHitGroupContribution);
};
