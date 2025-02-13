#pragma once

// Standard
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>

// Internal
#include "ArithmeticParameterOption.hpp"
#include "DefaultedParameterOption.hpp"
#include "DefaultedStringParameterOption.hpp"
#include "DirectoryParameterOption.hpp"
#include "FileParameterOption.hpp"
#include "GenomicOrientation.hpp"
#include "LogLevel.hpp"
#include "ParameterOption.hpp"

using namespace dataTypes;

struct GeneralOptions {
    static constexpr std::string_view optionsDescription =
        "RNAnue efficient data analysis for RNA–RNA interactomics.\nRun RNAnue with the subcall "
        "\"complete\" to execute all pipeline steps.\n\nMinimum call: RNAnue complete -t "
        "<treatment-dir> "
        "-o "
        "<output-dir> -f <feature-gff-file> --dbref <reference-genome-file>\nOr run RNAnue with a "
        "config "
        "file: RNAnue complete -c <config-file>\n\nGeneral Options"sv;

    static constexpr DirectoryParameterOption trtms{
        {.shortName = 't', .longName = "trtms"},
        "parent directory containing the raw read files of the treatments (required)"sv};

    static constexpr DirectoryParameterOption<true> ctrls{
        {.shortName = 's', .longName = "ctrls"},
        "parent directory containing the raw read files of the controls (optional)"sv};

    static constexpr ParameterOption<std::string> out{
        {.shortName = 'o', .longName = "outdir"},
        "output directory to which the results are saved (required)"sv};

    static constexpr DefaultedParameterOption<LogLevel, LogLevel::INFO> logLevel{
        {.shortName = std::nullopt, .longName = "loglevel"},
        "output directory to which the results are saved (required)"sv};

    static constexpr DefaultedParameterOption<size_t, 2> threads{
        {.shortName = 'p', .longName = "threads"}, "max number of threads to be used"sv};

    static constexpr FileParameterOption featuresPath{
        {.shortName = 'f', .longName = "features"},
        "annotation/features file in GFF/GTF format (required)"sv};

    static constexpr DefaultedStringParameterOption featureTypes{
        {.shortName = std::nullopt, .longName = "featuretypes"},
        "feature types to be considered for the analysis, can be specified as --featuretypes 'gene,rRNA' comma seperated values"sv,
        "transcript"sv};

    static constexpr DefaultedParameterOption<GenomicOrientation, GenomicOrientation::Value::BOTH>
        featureOrientation{
            {.shortName = std::nullopt, .longName = "orientation"},
            "orientation of the reads in relation to RNA sequences (strand-specific sequencing). Non strand-specific setting (both) disables nrgmax filtering. [same, opposite, both]"sv};

    static constexpr ArithmeticParameterOption<int, 100000,
                                               {.lowerBound = 1, .upperBound = 1000000}>
        chunkSize{{.shortName = std::nullopt, .longName = "chunksize"},
                  "number of reads processed per chunk in parallel"sv};

    static constexpr auto allOptions =
        std::make_tuple(trtms, ctrls, out, logLevel, threads, featuresPath, featureTypes,
                        featureOrientation, chunkSize);
};
