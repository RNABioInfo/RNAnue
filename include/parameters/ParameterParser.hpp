#pragma once

// Standard
#include <cstdlib>
#include <variant>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "AlignParameters.hpp"
#include "AnalyzeParameters.hpp"
#include "CompleteParameters.hpp"
#include "Config.hpp"
#include "DetectParameters.hpp"
#include "PostprocessParameters.hpp"
#include "PreprocessParameters.hpp"

namespace po = boost::program_options;

namespace pipelines {

class ParameterParser {
   public:
    using ParametersVariant =
        std::variant<CompleteParameters, preprocess::PreprocessParameters, align::AlignParameters,
                     detect::DetectParameters, analyze::AnalyzeParameters,
                     postprocess::PostprocessParameters>;

    static auto getParameters(int argc, const char* const argv[]) -> ParametersVariant;  // NOLINT

    ParameterParser() = delete;

   private:
    static auto parseParameters(int argc, const char* const argv[]) -> po::variables_map;  // NOLINT

    static void insertConfigFileParameters(po::variables_map& params);

    static auto getCommandLineOptions() -> po::options_description;
    static auto getConfigFileOptions() -> po::options_description;

    static auto getPositionalOptions() -> po::positional_options_description;

    static void printVersion();
};

}  // namespace pipelines
