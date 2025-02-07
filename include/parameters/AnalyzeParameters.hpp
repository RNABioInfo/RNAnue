#pragma once

// Standard
#include <cstddef>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "AnalyzeOptions.hpp"
#include "ClusteringParameters.hpp"
#include "GeneralParameters.hpp"
#include "GenomicStrandSpecificity.hpp"

namespace pipelines::analyze {

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
          clusterMergingStrandSpecificity(
              AnalyzeOptions::clusteringStrandSpecificity.extractValue(params)),
          maxClusterSelfOverlapFraction(AnalyzeOptions::maxSelfOverlap.extractValue(params)),
          padjThreshold(AnalyzeOptions::maxPadjValue.extractValue(params)),
          minimumClusterReadCount(AnalyzeOptions::minReadClusterCount.extractValue(params)) {};

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
        auto clusterOverlapFractionMin =
            AnalyzeOptions::clusterFractionOverlap.extractValue(params);

        if (clusterOverlapFractionMin) {
            return ShortestClusterOverlapFractionMergeParameter{clusterOverlapFractionMin.value()};
        }

        // Add one so that the distance threshold of 0 equals tolerance of one meaning blunt ended
        // clusters are merged by default.
        constexpr int DISTANCE_TO_TOLERANCE_OFFSET = 1;
        return ClusterOverlapToleranceMergeParameter{
            AnalyzeOptions::clusteringDistanceTolerance.extractValue(params) +
            DISTANCE_TO_TOLERANCE_OFFSET};
    }
};

}  // namespace pipelines::analyze
