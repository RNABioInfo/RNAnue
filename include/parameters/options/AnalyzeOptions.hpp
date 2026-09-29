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
    static constexpr std::string_view optionsDescription = "Analyze Pipeline";

    static constexpr DefaultedParameterOption<double, 0.1> maxSelfOverlap{
        {.shortName = std::nullopt, .longName = "max_self_overlap"},
        "maximum fractional overlap between two regions of a cluster"sv};

    static constexpr DefaultedParameterOption<GenomicStrandSpecificity,
                                              GenomicStrandSpecificity::UNSPECIFIC>
        clusteringStrandSpecificity{
            {.shortName = std::nullopt, .longName = "clustering_strand_specificity"},
            "whether to cluster strand specific or unspecific [specific, unspecific]"sv};

    static constexpr DefaultedParameterOption<int, 0> clusteringDistanceTolerance{
        {.shortName = std::nullopt, .longName = "clustering_distance"},
        "threshold distance at which two clusters are merged into a single combined cluster, default is to merge overlapping and blunt ended clusters (negative values can be used to implicitly expand search regions during clustering; mutually exclusive with --min_cluster_overlap)"sv};

    static constexpr ParameterOption<float, true> clusterFractionOverlap{
        {.shortName = std::nullopt, .longName = "min_cluster_overlap"},
        "minimal fractional overlap of the smaller cluster at which two clusters are merged into a single combined cluster, this option overwrites clustering_distance (mutually exclusive with --clustering_distance)"sv};

    static constexpr DefaultedParameterOption<double, 1.0> maxPadjValue{
        {.shortName = std::nullopt, .longName = "max_adjusted_p_value"},
        "adjusted abundance-corrected within-sample enrichment p-value threshold for outputting an interaction"sv};

    static constexpr DefaultedParameterOption<float, float{1.0}>
        minimumClusterTranscriptContribution{
            {.shortName = std::nullopt, .longName = "min_interaction_support"},
            "minimum weighted split-read contribution assigned to an interaction"sv};

    static constexpr ParameterOption<double, true> minimumSupportPerEffectiveBp{
        {.shortName = std::nullopt, .longName = "min_support_per_effective_bp"},
        "optional minimum weighted split-read support per effective coverage base pair"sv};

    static constexpr ParameterOption<size_t, true> maximumCoverageComponents{
        {.shortName = std::nullopt, .longName = "max_coverage_components"},
        "optional maximum number of thresholded coverage islands across both interaction arms"sv};

    static constexpr ParameterOption<double, true> minimumArmBalance{
        {.shortName = std::nullopt, .longName = "min_arm_balance"},
        "optional minimum ratio between weaker and stronger integrated arm coverage"sv};

    static constexpr auto allOptions =
        std::make_tuple(maxSelfOverlap, clusteringStrandSpecificity, clusteringDistanceTolerance,
                        clusterFractionOverlap, maxPadjValue, minimumClusterTranscriptContribution,
                        minimumSupportPerEffectiveBp, maximumCoverageComponents,
                        minimumArmBalance);
};
