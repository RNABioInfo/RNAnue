#pragma once

// Standard
#include <climits>
#include <cstddef>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <cstdint>

// Internal
#include "ClusteringParameters.hpp"
#include "Constants.hpp"
#include "GeneralParameters.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "ParameterValidator.hpp"

namespace pipelines::analyze {

using namespace constants::pipelines;
using namespace dataTypes;

class AnalyzeParameters : public GeneralParameters {
   public:
    ClusteringMergeParameterVariant clusteringParameters;
    GenomicStrandSpecificity clusterMergingStrandSpecificity;
    double maxClusterSelfOverlapFraction;
    double padjThreshold;
    size_t minimumClusterReadCount;

    AnalyzeParameters(const po::variables_map& params)
        : GeneralParameters(params),
          clusteringParameters(getClusteringMergeParameters(params)),
          clusterMergingStrandSpecificity(validateClusteringOrientation(params)),
          maxClusterSelfOverlapFraction(
              ParameterValidator::validateArithmetic(params, "maxselfoverlap", 0.0, 1.0)),
          padjThreshold(ParameterValidator::validateArithmetic(params, "padj", 0.0, 1.0)),
          minimumClusterReadCount(
              ParameterValidator::validateArithmetic<size_t>(params, "mincount", 1, SIZE_MAX)) {};

    static auto validateClusteringOrientation(const po::variables_map& params)
        -> GenomicStrandSpecificity {
        return params["clustmethod"].as<GenomicStrandSpecificity>();
    };

    auto getClusteringParameters() -> ClusteringParameters {
        return ClusteringParameters{
            .clusterMergeParameter = clusteringParameters,
            .clusterMergingStrandSpecificity = clusterMergingStrandSpecificity,
            .maxClusterSelfOverlapFraction = maxClusterSelfOverlapFraction,
            .minimumClusterReadCount = minimumClusterReadCount,
            .featureOrientation = featureOrientation};
    }

   private:
    static auto getClusteringMergeParameters(const po::variables_map& params)
        -> ClusteringMergeParameterVariant {
        if (params.contains("clustfrac")) {
            return ShortestClusterOverlapFractionMergeParameter{
                ParameterValidator::validateArithmetic<float>(params, "clustfrac",
                                                              clusterOverlapFractionMin, 1.0)};
        }

        // Add one so that the distance threshold of 0 equals tolerance of one meaning blunt ended
        // clusters are merged by default.
        constexpr int DISTANCE_TO_TOLERANCE_OFFSET = 1;
        return ClusterOverlapToleranceMergeParameter{
            ParameterValidator::validateArithmetic<int>(params, "clustdist", INT_MIN, INT_MAX) +
            DISTANCE_TO_TOLERANCE_OFFSET};
    }
};

}  // namespace pipelines::analyze
