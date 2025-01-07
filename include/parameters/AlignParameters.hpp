#pragma once

// Standard
#include <cstddef>
#include <cstdint>
#include <filesystem>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "GeneralParameters.hpp"
#include "ParameterValidator.hpp"

namespace po = boost::program_options;

namespace pipelines::align {

struct AlignParameters : public GeneralParameters {
    std::filesystem::path referenceGenome;
    bool multimapAlignments;
    size_t minLengthThreshold;
    size_t accuracy;
    size_t minimumFragmentScore;
    size_t minimumFragmentLength;
    size_t minimumSpliceCoverage;

    AlignParameters(const po::variables_map& params)
        : GeneralParameters(params),
          referenceGenome(ParameterValidator::validateFilePath(params, "dbref")),
          multimapAlignments(ParameterValidator::validateBool(params, "multimap")),
          minLengthThreshold(
              ParameterValidator::validateArithmetic<size_t>(params, "minlen", 0, 1000)),
          accuracy(ParameterValidator::validateArithmetic<size_t>(params, "accuracy", 0, 100)),
          minimumFragmentScore(
              ParameterValidator::validateArithmetic<size_t>(params, "minfragsco", 0, SIZE_MAX)),
          minimumFragmentLength(
              ParameterValidator::validateArithmetic<size_t>(params, "minfraglen", 0, SIZE_MAX)),
          minimumSpliceCoverage(
              ParameterValidator::validateArithmetic<size_t>(params, "minsplicecov", 0, 100)) {};
};

}  // namespace pipelines::align
