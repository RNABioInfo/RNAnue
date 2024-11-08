#pragma once

// Standard
#include <climits>
#include <cstddef>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "GeneralParameters.hpp"
#include "ParameterValidator.hpp"

namespace pipelines::analyze {

class AnalyzeParameters : public GeneralParameters {
   public:
    double maxOverlapFraction;
    int clusterDistanceThreshold;
    double padjThreshold;
    size_t minimumClusterReadCount;

    AnalyzeParameters(const po::variables_map& params)
        : GeneralParameters(params),
          maxOverlapFraction(
              ParameterValidator::validateArithmetic(params, "maxoverlap", 0.0, 1.0)),
          clusterDistanceThreshold(
              ParameterValidator::validateArithmetic(params, "clustdist", INT_MIN, INT_MAX)),
          padjThreshold(ParameterValidator::validateArithmetic(params, "padj", 0.0, 1.0)),
          minimumClusterReadCount(
              ParameterValidator::validateArithmetic<size_t>(params, "mincount", 1, SIZE_MAX)) {};
};

}  // namespace pipelines::analyze
