#pragma once

// Standard
#include <cstddef>
#include <optional>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "AnalyzeOptions.hpp"
#include "ClusteringParameters.hpp"
#include "GeneralParameters.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"

namespace pipelines::analyze {

using namespace dataTypes;

class AnalyzeParameters : public GeneralParameters {
   public:
    ClusteringMergeParameterVariant clusteringParameters;
    GenomicStrandSpecificity clusterMergingStrandSpecificity;
    double maxClusterSelfOverlapFraction;
    double padjThreshold;
    float minimumClusterTranscriptContribution;
    std::optional<double> minimumSupportPerEffectiveBp;
    std::optional<size_t> maximumCoverageComponents;
    std::optional<double> minimumArmBalance;

    AnalyzeParameters(const po::variables_map& params)
        : GeneralParameters(params),
          clusteringParameters(getClusteringMergeParameters(params)),
          clusterMergingStrandSpecificity(
              AnalyzeOptions::clusteringStrandSpecificity.extractValue(params)),
          maxClusterSelfOverlapFraction(AnalyzeOptions::maxSelfOverlap.extractValue(params)),
          padjThreshold(AnalyzeOptions::maxPadjValue.extractValue(params)),
          minimumClusterTranscriptContribution(
              AnalyzeOptions::minimumClusterTranscriptContribution.extractValue(params)),
          minimumSupportPerEffectiveBp(validateNonNegativeOptional(
              AnalyzeOptions::minimumSupportPerEffectiveBp.extractValue(params),
              AnalyzeOptions::minimumSupportPerEffectiveBp.getLongName())),
          maximumCoverageComponents(
              AnalyzeOptions::maximumCoverageComponents.extractValue(params)),
          minimumArmBalance(validateUnitIntervalOptional(
              AnalyzeOptions::minimumArmBalance.extractValue(params),
              AnalyzeOptions::minimumArmBalance.getLongName())) {};

    static auto validateClusteringOrientation(const po::variables_map& params)
        -> GenomicStrandSpecificity {
        return params["clustmethod"].as<GenomicStrandSpecificity>();
    };

    auto getClusteringParameters() -> ClusteringParameters {
        return ClusteringParameters{
            .clusterMergeParameter = clusteringParameters,
            .clusterMergingStrandSpecificity = clusterMergingStrandSpecificity,
            .maxClusterSelfOverlapFraction = maxClusterSelfOverlapFraction,
            .minimumClusterTrascriptContribution = minimumClusterTranscriptContribution,
            .featureOrientation = featureOrientation};
    }

   private:
    static auto getClusteringMergeParameters(const po::variables_map& params)
        -> ClusteringMergeParameterVariant {
        auto clusterOverlapFractionMin =
            AnalyzeOptions::clusterFractionOverlap.extractValue(params);

        if (clusterOverlapFractionMin) {
            return ShortestClusterOverlapFractionMergeParameter{clusterOverlapFractionMin.value()};
        }

        return ClusterOverlapToleranceMergeParameter{
            clusterDistanceToOverlapTolerance(
                AnalyzeOptions::clusteringDistanceTolerance.extractValue(params))};
    }

    [[nodiscard]] static auto validateNonNegativeOptional(std::optional<double> value,
                                                          const std::string& name)
        -> std::optional<double> {
        if (value && *value < 0.0) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(name,
                                                                " must be greater than or equal "
                                                                "to 0.");
        }

        return value;
    }

    [[nodiscard]] static auto validateUnitIntervalOptional(std::optional<double> value,
                                                           const std::string& name)
        -> std::optional<double> {
        if (value && (*value < 0.0 || *value > 1.0)) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(name,
                                                                " must be between 0 and 1.");
        }

        return value;
    }
};

}  // namespace pipelines::analyze
