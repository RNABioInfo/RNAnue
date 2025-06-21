#include "ParameterParser.hpp"

// Standard
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

// Boost
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "AlignOptions.hpp"
#include "AlignParameters.hpp"
#include "AnalyzeOptions.hpp"
#include "AnalyzeParameters.hpp"
#include "Closing.hpp"
#include "CompleteParameters.hpp"
#include "Constants.hpp"
#include "DetectOptions.hpp"
#include "DetectParameters.hpp"
#include "GeneralOptions.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "OtherOptions.hpp"
#include "ParameterOptions.hpp"
#include "PostprocessOptions.hpp"
#include "PreprocessOptions.hpp"
#include "PreprocessParameters.hpp"
#include "SubcallOptions.hpp"

namespace pipelines {

auto ParameterParser::getParameters(int argc, const char *const argv[])  // NOLINT
    -> ParameterParser::ParametersVariant {
    const auto params = parseParameters(argc, argv);

    const std::string subcall = params.at("subcall").as<std::string>();
    if (subcall == constants::pipelines::COMPLETE) {
        return CompleteParameters{params};
    }
    if (subcall == constants::pipelines::PREPROCESS) {
        return preprocess::PreprocessParameters{params};
    }
    if (subcall == constants::pipelines::ALIGN) {
        return align::AlignParameters{params};
    }
    if (subcall == constants::pipelines::DETECT) {
        return detect::DetectParameters{params};
    }
    if (subcall == constants::pipelines::ANALYZE) {
        return analyze::AnalyzeParameters{params};
    }
    // if (subcall == const)

    Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Unknown subcall: " + subcall);
    exit(EXIT_FAILURE);
}

auto ParameterParser::parseParameters(int argc,
                                      const char *const argv[]) -> po::variables_map {  // NOLINT
    const po::options_description commandLineOptions{getCommandLineOptions()};

    const po::positional_options_description positionalOptions{getPositionalOptions()};

    po::variables_map params;
    store(po::command_line_parser(argc, argv)
              .options(commandLineOptions)
              .positional(positionalOptions)
              .run(),
          params);

    notify(params);

    printVersion();

    if (params.at("version").as<bool>()) {
        Closing::printQuote();
        exit(EXIT_SUCCESS);
    }

    if (params.at("help").as<bool>()) {
        std::cout << commandLineOptions << "\n";
        Closing::printQuote();
        exit(EXIT_SUCCESS);
    }

    if (!params.contains("subcall")) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Please provide a subcall.");
    }

    Logger::setLogLevel(params.at("loglevel").as<LogLevel>());

    insertConfigFileParameters(params);

    return params;
}

void ParameterParser::insertConfigFileParameters(po::variables_map &params) {
    if (!params.contains("config")) {
        return;
    }

    const po::options_description configFileOptions{getConfigFileOptions()};

    const std::string configFilePath{params["config"].as<std::string>()};

    std::ifstream configIn{configFilePath};

    if (!configIn) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
            "Configuration file could not be opened!");
    }

    po::store(po::parse_config_file(configIn, configFileOptions), params);
    notify(params);
}

auto ParameterParser::getCommandLineOptions() -> po::options_description {
    const po::options_description generalOptions{ParameterOptions::getOptions<GeneralOptions>()};
    const po::options_description preprocessOptions{
        ParameterOptions::getOptions<PreprocessOptions>()};
    const po::options_description alignOptions{ParameterOptions::getOptions<AlignOptions>()};
    const po::options_description detectOptions{ParameterOptions::getOptions<DetectOptions>()};
    const po::options_description analyzeOptions{ParameterOptions::getOptions<AnalyzeOptions>()};
    const po::options_description postprocessOptions{
        ParameterOptions::getOptions<PostprocessOptions>()};
    const po::options_description otherOptions{ParameterOptions::getOptions<OtherOptions>()};
    const po::options_description subcallOptions{ParameterOptions::getOptions<SubcallOptions>()};

    po::options_description commandLineOptions{"Command line options"};

    commandLineOptions.add(generalOptions)
        .add(preprocessOptions)
        .add(alignOptions)
        .add(detectOptions)
        .add(analyzeOptions)
        .add(postprocessOptions)
        .add(otherOptions)
        .add(subcallOptions);

    return commandLineOptions;
}

auto ParameterParser::getConfigFileOptions() -> po::options_description {
    const po::options_description generalOptions{ParameterOptions::getOptions<GeneralOptions>()};
    const po::options_description preprocessOptions{
        ParameterOptions::getOptions<PreprocessOptions>()};
    const po::options_description alignOptions{ParameterOptions::getOptions<AlignOptions>()};
    const po::options_description detectOptions{ParameterOptions::getOptions<DetectOptions>()};
    const po::options_description analyzeOptions{ParameterOptions::getOptions<AnalyzeOptions>()};

    po::options_description configFileOptions{"Config file options"};

    configFileOptions.add(generalOptions)
        .add(preprocessOptions)
        .add(alignOptions)
        .add(detectOptions)
        .add(analyzeOptions);

    return configFileOptions;
}

auto ParameterParser::getPositionalOptions() -> po::positional_options_description {
    po::positional_options_description positionalOptions;
    positionalOptions.add("subcall", 1);

    return positionalOptions;
}

void ParameterParser::printVersion() {
    const std::string versionString =
        "RNAnue v" + std::to_string(RNAnue_VERSION_MAJOR) + "." +
        std::to_string(RNAnue_VERSION_MINOR) + "." + std::to_string(RNAnue_VERSION_PATCH) + " - " +
        "Detect RNA-RNA interactions from Direct-Duplex-Detection (DDD) data.";

    Logger::log(versionString);
}

}  // namespace pipelines
