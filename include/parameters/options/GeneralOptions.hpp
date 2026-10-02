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
        "\"complete\" to execute all pipeline steps.\n\nMinimum call: RNAnue complete -T "
        "<treatment-dir> "
        "-o "
        "<output-dir> -f <feature-gff-file> --reference_genome <reference-genome-file>\nOr run "
        "RNAnue with a "
        "config "
        "file: RNAnue complete -c <config-file>\n\nGeneral Options"sv;

    static constexpr DirectoryParameterOption trtms{
        {.shortName = 'T', .longName = "treatment_dir"},
        "parent directory containing the raw read files of the treatments (required)"sv};

    static constexpr DirectoryParameterOption<true> ctrls{
        {.shortName = 'C', .longName = "control_dir"},
        "parent directory containing the raw read files of the controls (optional)"sv};

    static constexpr ParameterOption<std::string> out{
        {.shortName = 'o', .longName = "output_dir"},
        "output directory to which the results are saved (required)"sv};

    static constexpr DefaultedParameterOption<LogLevel, LogLevel::INFO> logLevel{
        {.shortName = std::nullopt, .longName = "log_level"},
        "logging verbosity [debug, info, warning, error]"sv};

    static constexpr ArithmeticParameterOption<size_t, 2, {.lowerBound = 2, .upperBound = 1000}>
        threads{{.shortName = 't', .longName = "threads"}, "max number of threads to be used"sv};

    static constexpr FileParameterOption featuresPath{
        {.shortName = 'f', .longName = "features"},
        "annotation/features file in GFF/GTF format (required)"sv};

    static constexpr DefaultedStringParameterOption featureTypes{
        {.shortName = std::nullopt, .longName = "feature_types"},
        "feature types to be considered for the analysis, can be specified as --feature_types 'gene,rRNA' comma seperated values"sv,
        "transcript"sv};

    static constexpr DefaultedParameterOption<GenomicOrientation, GenomicOrientation::Value::BOTH>
        featureOrientation{
            {.shortName = std::nullopt, .longName = "orientation"},
            "orientation of the reads in relation to RNA sequences (strand-specific sequencing). Non strand-specific setting (both) disables max_hybridization_energy filtering. [same, opposite, both]"sv};
    static constexpr DefaultedParameterOption<bool, false> maskMultiCopyGenes{
        {.shortName = std::nullopt,
         .longName = "mask_multicopy_genes"},
        "mask multi-copy genes, only keeping one transcript per gene"sv};

    static constexpr ArithmeticParameterOption<double, 0.99, {.lowerBound = 0.0, .upperBound = 1.0}>
        minMultiCopyIdentity{{.shortName = std::nullopt, .longName = "min_copy_identity"},
                             "minimum identity of two transcripts to be considered a copy"sv};

    static constexpr ArithmeticParameterOption<int, 100000,
                                               {.lowerBound = 1, .upperBound = 1000000}>
        chunkSize{{.shortName = std::nullopt, .longName = "chunk_size"},
                  "number of reads processed per chunk in parallel"sv};

    static constexpr auto allOptions =
        std::make_tuple(trtms, ctrls, out, logLevel, threads, featuresPath, featureTypes,
                        featureOrientation, maskMultiCopyGenes, minMultiCopyIdentity, chunkSize);
};
