#pragma once

// Standard
#include <algorithm>
#include <cctype>
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
#include "GeneralOptions.hpp"
#include "GenomicOrientation.hpp"
#include "LogLevel.hpp"

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

    // TODO: Write standalone featureTypes from string
    GeneralParameters(const po::variables_map& params)
        : treatmentsDir(GeneralOptions::trtms.extractValue(params)),
          controlDir(GeneralOptions::ctrls.extractValue(params)),
          outputDir(GeneralOptions::out.extractValue(params)),
          featuresInPath(GeneralOptions::featuresPath.extractValue(params)),
          featureTypes(validateFeatureTypes(params)),
          featureOrientation(GeneralOptions::featureOrientation.extractValue(params)),
          logLevel(GeneralOptions::logLevel.extractValue(params)),
          threadCount(GeneralOptions::threads.extractValue(params)),
          chunkSize(GeneralOptions::chunkSize.extractValue(params)) {};

   private:
    static auto validateFeatureTypes(const po::variables_map& params)
        -> std::unordered_set<std::string> {
        const auto featureTypesString = params["featuretypes"].as<std::string>();

        std::unordered_set<std::string> uniqueIncludedFeatures;

        std::stringstream stringStream(featureTypesString);
        std::string str;
        while (getline(stringStream, str, ',')) {
            str.erase(remove_if(str.begin(), str.end(), isspace), str.end());  // NOLINT
            uniqueIncludedFeatures.insert(str);
        }

        return uniqueIncludedFeatures;
    }
};
