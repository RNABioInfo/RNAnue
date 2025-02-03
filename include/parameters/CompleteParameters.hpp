#pragma once

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "AlignParameters.hpp"
#include "AnalyzeParameters.hpp"
#include "DetectParameters.hpp"
#include "PreprocessParameters.hpp"

namespace po = boost::program_options;

namespace pipelines {

struct CompleteParameters {
    preprocess::PreprocessParameters preprocessParameters;
    align::AlignParameters alignParameters;
    detect::DetectParameters detectParameters;
    analyze::AnalyzeParameters analyzeParameters;

    CompleteParameters(const po::variables_map &params)
        : preprocessParameters(params),
          alignParameters(params),
          detectParameters(params),
          analyzeParameters(params) {};
};

}  // namespace pipelines
