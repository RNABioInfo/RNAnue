#pragma once

// Standard
#include <sys/stat.h>

#include <cstddef>
#include <climits>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>

// Internal
#include "ArithmeticParameterOption.hpp"
#include "DefaultedParameterOption.hpp"
#include "DefaultedStringParameterOption.hpp"
#include "ParameterOption.hpp"

struct AlignOptions {
    static constexpr std::string_view optionsDescription{"Align Pipeline"sv};

    static constexpr ParameterOption<std::string> refGenome{
        {.shortName = 'r', .longName = "reference_genome"}, "reference genome (.fasta) (required)"sv};

    static constexpr DefaultedStringParameterOption aligner{
        {.shortName = 'a', .longName = "aligner"},
        "alignment backend [segemehl, star]; STAR execution is not implemented yet"sv,
        "segemehl"sv};

    static constexpr DefaultedParameterOption<bool, true> allowMultimap{
        {.shortName = std::nullopt,
         .longName = "allow_multimapping",
         .inverseLongName = "no_multimapping",
         .inverseDescription = "disable multimapping alignments"},
        "consider multimapping alignments."sv};

    static constexpr ArithmeticParameterOption<size_t, 90, {.lowerBound = 0, .upperBound = 100}>
        accuracy{{.shortName = std::nullopt, .longName = "seg_alignment_accuracy"},
                 "segemehl-specific minimum percentage of read matches"sv};

    static constexpr ArithmeticParameterOption<size_t, 18, {.lowerBound = 0, .upperBound = 100}>
        minFragmentScore{{.shortName = std::nullopt, .longName = "seg_min_fragment_score"},
                         "segemehl-specific minimum score of a spliced fragment"sv};

    static constexpr ArithmeticParameterOption<size_t, 20,
                                               {.lowerBound = 0, .upperBound = SIZE_MAX}>
        minAlignLength{{.shortName = std::nullopt, .longName = "min_alignment_length"},
                       "minimum total length of the aligned fraction"sv};

    static constexpr ArithmeticParameterOption<size_t, 10,
                                               {.lowerBound = 0, .upperBound = SIZE_MAX}>
        minFragmentLength{{.shortName = std::nullopt, .longName = "min_fragment_length"},
                          "minimum length of a spliced fragment"sv};

    static constexpr ArithmeticParameterOption<size_t, 80, {.lowerBound = 0, .upperBound = 100}>
        minSpliceCoverage{{.shortName = std::nullopt, .longName = "min_split_coverage"},
                          "minimum coverage for spliced transcripts"sv};

    // Parameter preparation only: these controls do not enable STAR execution.
    static constexpr ArithmeticParameterOption<int, 10, {.lowerBound = 1, .upperBound = INT_MAX}>
        starMaxMultimaps{{.shortName = std::nullopt, .longName = "star_max_multimaps"},
                        "future STAR outFilterMultimapNmax/chimMultimapNmax; effective cap is 1 when multimapping is disabled"sv};

    static constexpr ParameterOption<int, true> starMinJunctionOverhang{
        {.shortName = std::nullopt, .longName = "star_min_junction_overhang"},
        "future STAR chimJunctionOverhangMin (>=1); omitted: inherit min_fragment_length"sv};

    static constexpr ArithmeticParameterOption<int, 3, {.lowerBound = 0, .upperBound = INT_MAX}>
        starMaxSegmentGap{{.shortName = std::nullopt, .longName = "star_max_segment_gap"},
                          "future STAR chimSegmentReadGapMax"sv};

    static constexpr ArithmeticParameterOption<int, 10, {.lowerBound = 0, .upperBound = INT_MAX}>
        starMinNonchimericScoreDrop{
            {.shortName = std::nullopt, .longName = "star_min_nonchimeric_score_drop"},
            "future STAR chimNonchimScoreDropMin"sv};

    static constexpr ArithmeticParameterOption<int, 30, {.lowerBound = 0, .upperBound = INT_MAX}>
        starMaxChimericScoreDrop{
            {.shortName = std::nullopt, .longName = "star_max_chimeric_score_drop"},
            "future STAR chimScoreDropMax"sv};

    static constexpr ArithmeticParameterOption<int, 0, {.lowerBound = 0, .upperBound = INT_MAX}>
        starMaxIntronLength{{.shortName = std::nullopt, .longName = "star_max_intron_length"},
                           "future STAR alignIntronMax; 0 uses STAR's automatic bound"sv};

    static constexpr auto starOptions =
        std::make_tuple(starMaxMultimaps, starMinJunctionOverhang, starMaxSegmentGap,
                        starMinNonchimericScoreDrop, starMaxChimericScoreDrop, starMaxIntronLength);

    static constexpr auto allOptions = std::tuple_cat(
        std::make_tuple(refGenome, aligner, allowMultimap, accuracy, minFragmentScore,
                        minAlignLength, minFragmentLength, minSpliceCoverage), starOptions);
};
