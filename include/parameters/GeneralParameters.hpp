#pragma once

// Standard
#include <algorithm>
#include <climits>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_set>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "GenomicOrientation.hpp"
#include "Logger.hpp"
#include "ParameterValidator.hpp"

namespace po = boost::program_options;
class GeneralParameters {
   public:
    std::filesystem::path treatmentsDir;
    std::optional<std::filesystem::path> controlDir;
    std::filesystem::path outputDir;

    std::filesystem::path featuresInPath;
    std::unordered_set<std::string> featureTypes;
    dataTypes::GenomicOrientation featureOrientation;

    LogLevel logLevel;

    size_t threadCount;
    size_t chunkSize;

    GeneralParameters(const po::variables_map& params)
        : treatmentsDir(ParameterValidator::validateDirectory(params, "trtms")),
          controlDir(validateControlDir(params)),
          outputDir(ParameterValidator::validateDirectory(params, "outdir")),
          featuresInPath(ParameterValidator::validateFilePath(params, "features")),
          featureTypes(validateFeatureTypes(params)),
          featureOrientation(validateFeatureOrientation(params)),
          logLevel(validateLogLevel(params)),
          threadCount(ParameterValidator::validateArithmetic(params, "threads", 1, INT_MAX)),
          chunkSize(ParameterValidator::validateArithmetic(params, "chunksize", 1, INT_MAX)) {};

   private:
    static auto validateControlDir(const po::variables_map& params)
        -> std::optional<std::filesystem::path> {
        if ((params.count("ctrls") != 0U) && !params["ctrls"].as<std::string>().empty()) {
            return ParameterValidator::validateDirectory(params, "ctrls");
        }

        return std::nullopt;
    }

    static auto validateFeatureTypes(const po::variables_map& params)
        -> std::unordered_set<std::string> {
        const auto featureTypesString = params["featuretypes"].as<std::string>();

        std::unordered_set<std::string> uniqueIncludedFeatures;

        std::stringstream stringStream(featureTypesString);
        std::string str;
        while (getline(stringStream, str, ',')) {
            str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
            uniqueIncludedFeatures.insert(str);
        }

        return uniqueIncludedFeatures;
    }

    static auto validateFeatureOrientation(const po::variables_map& params)
        -> dataTypes::GenomicOrientation {
        return params["orientation"].as<dataTypes::GenomicOrientation>();
    }

    static auto validateLogLevel(const po::variables_map& params) -> LogLevel {
        const std::string logLevelStr = params["loglevel"].as<std::string>();

        if (logLevelStr == "debug" || logLevelStr == "DEBUG") {
            return LogLevel::DEBUG;
        }
        if (logLevelStr == "info" || logLevelStr == "INFO") {
            return LogLevel::INFO;
        }
        if (logLevelStr == "warning" || logLevelStr == "WARNING") {
            return LogLevel::WARNING;
        }
        if (logLevelStr == "error" || logLevelStr == "ERROR") {
            return LogLevel::ERROR;
        }
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Invalid log level specified.");
        return LogLevel::INFO;
    }
};
