#pragma once

// Standard
#include <cstddef>
#include <optional>
#include <string_view>
#include <tuple>

// Boost
#include <boost/program_options/options_description.hpp>

// Internal
#include "ArithmeticParameterOption.hpp"
#include "DefaultedParameterOption.hpp"
#include "DefaultedStringParameterOption.hpp"
#include "ValueBounds.hpp"

struct PreprocessOptions {
    static constexpr std::string_view optionsDescription = "Preprocess Pipeline"sv;

    // Boolean parameters (using DefaultedParameterOption)
    static constexpr DefaultedParameterOption<bool, true> enablePreprocess{
        {.shortName = std::nullopt, .longName = "preprocess"},
        "whether to include preprocessing of the raw reads in the workflow of RNAnue"sv};

    static constexpr DefaultedParameterOption<bool, true> enableDeduplicate{
        {.shortName = std::nullopt, .longName = "deduplicate"},
        "whether to remove duplicate reads based on the sequence"sv};

    static constexpr DefaultedParameterOption<bool, true> trimPolyG{
        {.shortName = std::nullopt, .longName = "trimpolyg"},
        "whether to trim high quality polyG tails from the reads. Applicable for Illumina NextSeq reads"sv};

    static constexpr ArithmeticParameterOption<
        std::size_t, 5, ValueBounds<std::size_t>{.lowerBound = 0, .upperBound = 100}>
        minPolyGCount{{.shortName = std::nullopt, .longName = "minpolygcount"},
                      "minimum count of consecutive G's to be cut off by trimpolyg"sv};

    // String parameters for adapter sequences.
    static constexpr DefaultedStringParameterOption adpt5f{
        {.shortName = std::nullopt, .longName = "adpt5f"},
        "single sequence or file [.fasta] of the adapter sequences to be removed from the 5' end (forward read if PE)"sv,
        {""sv}};

    static constexpr DefaultedStringParameterOption adpt5r{
        {.shortName = std::nullopt, .longName = "adpt5r"},
        "single sequence or file [.fasta] of the adapter sequences to be removed from the 5' end of the reverse read (PE only)"sv,
        ""sv};

    static constexpr DefaultedStringParameterOption adpt3f{
        {.shortName = std::nullopt, .longName = "adpt3f"},
        "single sequence or file [.fasta] of the adapter sequences to be removed from the 3' end (forward read if PE)"sv,
        ""sv};

    static constexpr DefaultedStringParameterOption adpt3r{
        {.shortName = std::nullopt, .longName = "adpt3r"},
        "single sequence or file [.fasta] of the adapter sequences to be removed from the 3' end of the reverse read (PE only)"sv,
        ""sv};

    // Double (arithmetic) parameters.
    // For mtrim: a mismatch rate (between 0.0 and 1.0).
    static constexpr ArithmeticParameterOption<
        double, 0.05, ValueBounds<double>{.lowerBound = 0.0, .upperBound = 1.0}>
        mtrim{{.shortName = std::nullopt, .longName = "mtrim"},
              "rate of mismatches allowed when aligning adapters to sequences"sv};

    // For minovltrim: minimum adapter-read overlap.
    static constexpr ArithmeticParameterOption<
        std::size_t, 5, ValueBounds<std::size_t>{.lowerBound = 0, .upperBound = 1000}>
        minOvlTrim{{.shortName = std::nullopt, .longName = "minovltrim"},
                   "minimum length of overlap between adapter and read"sv};

    // For minqual: the lower limit for the mean quality.
    static constexpr ArithmeticParameterOption<
        std::size_t, 20, ValueBounds<std::size_t>{.lowerBound = 0, .upperBound = 100}>
        minQual{{.shortName = 'q', .longName = "minqual"},
                "lower limit for the mean quality (Phred Quality Score) of the reads"sv};

    // For minlen: the minimum read length.
    static constexpr ArithmeticParameterOption<
        std::size_t, 15, ValueBounds<std::size_t>{.lowerBound = 0, .upperBound = 1000}>
        minLen{{.shortName = 'l', .longName = "minlen"},
               "minimum length of the reads (default: 15)"sv};

    // For wqual: the minimum window quality.
    static constexpr ArithmeticParameterOption<
        std::size_t, 20, ValueBounds<std::size_t>{.lowerBound = 0, .upperBound = 100}>
        wqual{{.shortName = std::nullopt, .longName = "wqual"},
              "minimum mean quality for each window (Phred Quality Score)"sv};

    // For wtrim: the window size for quality trimming.
    static constexpr ArithmeticParameterOption<
        std::size_t, 0, ValueBounds<std::size_t>{.lowerBound = 0, .upperBound = 1000}>
        wtrim{
            {.shortName = std::nullopt, .longName = "wtrim"},
            "window size for quality trimming from 3' end. Selecting '0' will not apply quality trimming"sv};

    // For minovl: the minimal overlap for merging paired-end reads.
    static constexpr ArithmeticParameterOption<
        std::size_t, 5, ValueBounds<std::size_t>{.lowerBound = 0, .upperBound = 1000}>
        minOvl{{.shortName = std::nullopt, .longName = "minovl"},
               "minimal overlap to merge paired-end reads"sv};

    // For mmerge: the mismatch rate when merging paired-end reads.
    static constexpr ArithmeticParameterOption<
        double, 0.05, ValueBounds<double>{.lowerBound = 0.0, .upperBound = 1.0}>
        mmerge{{.shortName = std::nullopt, .longName = "mmerge"},
               "rate of mismatches allowed when merging paired end reads"sv};

    static constexpr auto allOptions = std::make_tuple(
        enablePreprocess, enableDeduplicate, trimPolyG, minPolyGCount, adpt5f, adpt5r, adpt3f,
        adpt3r, mtrim, minOvlTrim, minQual, minLen, wqual, wtrim, minOvl, mmerge);
};
