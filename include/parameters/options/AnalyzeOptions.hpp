#pragma once

// Standard
#include <cstddef>
#include <optional>
#include <string_view>
#include <tuple>

// Internal
#include "DefaultedParameterOption.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "ParameterOption.hpp"

using namespace dataTypes;

struct AnalyzeOptions {
    static constexpr std::string_view optionsDescription = "Analyze Pipeline"sv;

    static constexpr DefaultedParameterOption<double, 0.1> maxSelfOverlap{
        {.shortName = std::nullopt, .longName = "maxselfoverlap"},
        "maximum fractional overlap between two regions of a cluster"sv};

    static constexpr DefaultedParameterOption<GenomicStrandSpecificity,
                                              GenomicStrandSpecificity::UNSPECIFIC>
        clusteringStrandSpecificity{
            {.shortName = std::nullopt, .longName = "clustmethod"},
            "whether to cluster strand specific or unspecific [specific, unspecific]"sv};

    static constexpr DefaultedParameterOption<int, 0> clusteringDistanceTolerance{
        {.shortName = std::nullopt, .longName = "clustdist"},
        "threshold distance at which two clusters are merged into a single combined cluster, default is to only merge overlapping and blunt ended clusters (mutually exclusive with --clustfrac)"sv};

    static constexpr ParameterOption<float, true> clusterFractionOverlap{
        {.shortName = std::nullopt, .longName = "clustfrac"},
        "minimal fractional overlap of the smaller cluster at which two clusters are merged into a single combined cluster, this option overwrites clustdist (mutually exclusive with --clustdist)"sv};

    static constexpr DefaultedParameterOption<double, 1.0> maxPadjValue{
        {.shortName = std::nullopt, .longName = "padj"},
        "adjusted p-value threshold for outputting an interaction"sv};

    static constexpr DefaultedParameterOption<float, float{1.0}>
        minimumClusterTranscriptContribution{
            {.shortName = std::nullopt, .longName = "mincount"},
            "minimum number of scored transcripts assigned to an interaction"sv};

    static constexpr auto allOptions =
        std::make_tuple(maxSelfOverlap, clusteringStrandSpecificity, clusteringDistanceTolerance,
                        clusterFractionOverlap, maxPadjValue, minimumClusterTranscriptContribution);
};
