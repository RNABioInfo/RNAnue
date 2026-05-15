#pragma once

// Standard
#include <sys/stat.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>

// Internal
#include "ArithmeticParameterOption.hpp"
#include "DefaultedParameterOption.hpp"
#include "ParameterOption.hpp"

struct AlignOptions {
    static constexpr std::string_view optionsDescription{"Align Pipeline"sv};

    static constexpr ParameterOption<std::string> refGenome{
        {.shortName = std::nullopt, .longName = "dbref"}, "reference genome (.fasta) (required)"sv};

    static constexpr DefaultedParameterOption<bool, true> allowMultimap{
        {.shortName = std::nullopt,
         .longName = "multimap",
         .inverseLongName = "no-multimap",
         .inverseDescription = "disable multimapping alignments"},
        "consider multimapping alignments."sv};

    static constexpr ArithmeticParameterOption<size_t, 90, {.lowerBound = 0, .upperBound = 100}>
        accuracy{{.shortName = std::nullopt, .longName = "accuracy"},
                 "minimum percentage of read matches"sv};

    static constexpr ArithmeticParameterOption<size_t, 18, {.lowerBound = 0, .upperBound = 100}>
        minFragmentScore{{.shortName = std::nullopt, .longName = "minfragsco"},
                         "minimum score of a spliced fragment"sv};

    static constexpr ArithmeticParameterOption<size_t, 20,
                                               {.lowerBound = 0, .upperBound = SIZE_MAX}>
        minAlignLength{{.shortName = std::nullopt, .longName = "minalignlen"},
                       "minimum total length of the aligned fraction"sv};

    static constexpr ArithmeticParameterOption<size_t, 10,
                                               {.lowerBound = 0, .upperBound = SIZE_MAX}>
        minFragmentLength{{.shortName = std::nullopt, .longName = "minfraglen"},
                          "minimum length of a spliced fragment"sv};

    static constexpr ArithmeticParameterOption<size_t, 80, {.lowerBound = 0, .upperBound = 100}>
        minSpliceCoverage{{.shortName = std::nullopt, .longName = "minsplicecov"},
                          "minimum coverage for spliced transcripts"sv};

    static constexpr auto allOptions =
        std::make_tuple(refGenome, allowMultimap, accuracy, minFragmentScore, minAlignLength,
                        minFragmentLength, minSpliceCoverage);
};
