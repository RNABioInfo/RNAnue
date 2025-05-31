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
        maxPrimaryAlignmentCount{{.shortName = std::nullopt, .longName = "maxprim"},
                                 "maximum count of primary alignments per read"sv};

    static constexpr ArithmeticParameterOption<size_t, 0, {.lowerBound = 0, .upperBound = 100}>
        minMappingQuality{{.shortName = std::nullopt, .longName = "mapqmin"},
                          "minimum mapping quality for reads to be considered in the analysis"sv};

    static constexpr ArithmeticParameterOption<double, 0.5, {.lowerBound = 0.0, .upperBound = 1.0}>
        minComplementarity{{.shortName = std::nullopt, .longName = "cmplmin"},
                           "complementarity cutoff for split reads"sv};

    static constexpr ArithmeticParameterOption<double, 0.1, {.lowerBound = 0.0, .upperBound = 1.0}>
        siteLengthRatio{{.shortName = std::nullopt, .longName = "sitelenratio"},
                        "fraction of shortest record match cutoff"sv};

    static constexpr ArithmeticParameterOption<size_t, 20,
                                               {.lowerBound = 0, .upperBound = SIZE_MAX}>
        minDetectLength{{.shortName = std::nullopt, .longName = "mindetectlen"},
                        "minimum fragment length after clipping"sv};

    static constexpr DefaultedParameterOption<double, 0.0> maxEnergy{
        {.shortName = std::nullopt, .longName = "nrgmax"},
        "hybridization energy cutoff for split reads"sv};

    static constexpr DefaultedParameterOption<bool, false> excludeSoftClipping{
        {.shortName = std::nullopt, .longName = "exclclipping"},
        "exclude soft clipping from the alignments"sv};

    static constexpr DefaultedParameterOption<bool, false> filterSplicing{
        {.shortName = std::nullopt, .longName = "splicing"},
        "splicing events are removed in the detection of split reads"sv};

    static constexpr DefaultedParameterOption<bool, true> allowAltSplicing{
        {.shortName = std::nullopt, .longName = "altsplice"},
        "remove alternative splicing events"sv};

    static constexpr DefaultedParameterOption<int, 5> splicingTolerance{
        {.shortName = std::nullopt, .longName = "spltol"}, "tolerance for splicing events"sv};

    static constexpr DefaultedParameterOption<bool, false> includeWobble{
        {.shortName = std::nullopt, .longName = "includewobble"},
        "[EXPERIMENTAL] wobble base pairs are allowed in crosslinking site evaluation"sv};

    static constexpr ArithmeticParameterOption<double, 0.1, {.lowerBound = 0, .upperBound = 1}>
        minHitGroupContribution{
            {.shortName = std::nullopt, .longName = "mincontr"},
            "minimum contribution of alignment within all alignments of this read"sv};

    static constexpr auto allOptions = std::make_tuple(
        maxPrimaryAlignmentCount, minMappingQuality, minComplementarity, siteLengthRatio,
        minDetectLength, maxEnergy, excludeSoftClipping, filterSplicing, allowAltSplicing,
        splicingTolerance, includeWobble, minHitGroupContribution);
};
