#include "FeatureParser.hpp"

// Standard
#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Internal
#include "Constants.hpp"
#include "FileType.hpp"
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"

using namespace constants::annotation;

namespace annotation {

/** @brief Constructs a FeatureParser object.
 *
 * This constructor initializes a FeatureParser object with the specified included features and
 * feature ID flag.
 *
 * @param includedFeatures The set of included features. Default is exon.
 *        Possible values include exon, gene, lncRNA, and other.
 * @param featureIDFlag The optional feature ID flag. Default is empty.
 *       Default is ID for GFF files and gene_id for GTF files.
 */
FeatureParser::FeatureParser(const std::unordered_set<std::string> &includedFeatures,
                             const std::optional<std::string> &featureIDFlag)
    : includedFeatures(includedFeatures), featureIDFlag(featureIDFlag) {}

auto FeatureParser::parse(const fs::path &featureFilePath,
                          const ReferenceIDToIndexMap &referenceToIndex) const -> FeatureMap {
    FileType fileType = getFileType(featureFilePath);

    return iterateFeatureFile(featureFilePath, fileType, referenceToIndex);
}

auto FeatureParser::getFileType(const fs::path &featureFilePath) -> FileType {
    std::ifstream file(featureFilePath.string());
    if (!file) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not open file: " +
                                                            featureFilePath.string());
    }

    std::string line;
    std::getline(file, line);

    if (line.starts_with("##gff-version")) {
        return FileType(FileType::GFF);
    }
    if (line.starts_with("##gtf-version")) {
        return FileType(FileType::GTF);
    }

    throw std::runtime_error(
        "Annotation file type not supported. First line does not contain "
        "##gff-version or ##gtf-version");
}

auto FeatureParser::iterateFeatureFile(const fs::path &featureFilePath, const FileType fileType,
                                       const ReferenceIDToIndexMap &referenceToIndex) const
    -> FeatureMap {
    FeatureMap featureMap;

    std::ifstream file(featureFilePath.string());

    if (!file.is_open()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not open file: ",
                                                            featureFilePath.string());
    }

    size_t parsedFeatures = 0;
    std::unordered_set<std::string> featureGroups;

    for (std::string line; std::getline(file, line);) {
        if (line[0] == '#') {
            continue;
        }

        const auto tokens = getTokens(line, includedFeatures);

        if (!isValidFeature(tokens)) {
            continue;
        }

        const auto &tokens_v = tokens.value();

        const std::string &referenceID = tokens_v[0];

        if (!referenceToIndex.contains(referenceID)) {
            Logger::log<LogLevel::WARNING>(
                "Feature reference id not found: ", referenceID,
                ", ensure the annotion and reference genome use the same reference IDs.");

            continue;
        }

        const int referenceIDIndex = referenceToIndex.at(referenceID);
        const std::string &featureType = tokens_v[2];

        // Convert to zero-based half-open index
        int startPosition = std::stoi(tokens_v[3]) - 1;
        int endPosition = std::stoi(tokens_v[4]);

        const auto attributes = getAttributes(fileType, tokens_v[8]);

        auto getAttribute = [&attributes](const std::string &key) -> std::optional<std::string> {
            const auto iterator = attributes.find(key);
            if (iterator == attributes.end()) {
                return std::nullopt;
            }
            return iterator->second;
        };

        const std::string featureIDFlag = this->featureIDFlag.value_or(fileType.defaultIDKey());
        const auto identifier = getAttribute(featureIDFlag);

        if (!identifier.has_value()) {
            Logger::log<LogLevel::WARNING>("Could not find identifier in GFF file");
            continue;
        }

        const std::optional<std::string> geneName =
            getAttribute(annotation::FileType::defaultGeneNameKey());

        GenomicRegion genomicRegion{referenceIDIndex,
                                    {.startPosition = startPosition, .endPosition = endPosition},
                                    getGenomicStrand(tokens_v[strandTokenColumn][0])};

        featureMap[referenceIDIndex].emplace_back(featureType, genomicRegion, identifier.value(),
                                                  getAttribute(fileType.defaultGroupKey()),
                                                  geneName);

        ++parsedFeatures;

        if (featureMap[referenceIDIndex].back().getGroupID().has_value()) {
            featureGroups.insert(featureMap[referenceIDIndex].back().getGroupID().value());
        }
    }

    const std::string includedFeatureTypes =
        std::accumulate(includedFeatures.begin(), includedFeatures.end(), std::string(),
                        [](const std::string &lhs, const std::string &rhs) {
                            return lhs.empty() ? rhs : lhs + ", " + rhs;
                        });

    if (parsedFeatures != 0) {
        const std::string featureGroupLog =
            featureGroups.empty()
                ? ""
                : " Found " + std::to_string(featureGroups.size()) + " feature groups.";

        Logger::log("Parsed ", std::to_string(parsedFeatures),
                    " features of type: ", includedFeatureTypes, ".", featureGroupLog);
    } else {
        Logger::log<LogLevel::WARNING>(
            "No features parsed. Check your feature file and included features flag.");
    }

    return featureMap;
}

auto FeatureParser::getTokens(const std::string &line,
                              const std::unordered_set<std::string> &includedFeatures)
    -> std::optional<std::vector<std::string>> {
    std::vector<std::string> tokens;
    std::istringstream issLine(line);

    for (std::string token; std::getline(issLine, token, '\t');) {
        // Checks at the third token if the feature is included
        if (tokens.size() == 2 && !includedFeatures.contains(token)) {
            return std::nullopt;
        }

        tokens.push_back(token);
    }

    return tokens;
}

auto FeatureParser::getAttributes(const annotation::FileType fileType,
                                  const std::string &attributes)
    -> std::unordered_map<std::string, std::string> {
    std::unordered_map<std::string, std::string> attributeMap;
    std::istringstream issAttr(attributes);

    const char delim = annotation::FileType::attributeDelimiter();

    for (std::string attribute; std::getline(issAttr, attribute, delim);) {
        const auto keyPosition = attribute.find(fileType.attributeAssignment());

        if (keyPosition == std::string::npos) {
            continue;
        }

        std::string key = attribute.substr(0, keyPosition);
        std::string value = attribute.substr(keyPosition + 1);

        if (key.empty() || value.empty()) {
            continue;
        }

        if (fileType == FileType::GTF) {
            const auto keyRet = std::ranges::remove(key, ' ');
            key.erase(keyRet.begin(), keyRet.end());
            const auto attributeRet = std::ranges::remove(value, '\"');
            attribute.erase(attributeRet.begin(), attributeRet.end());
        }

        attributeMap[key] = value;
    }

    return attributeMap;
}

constexpr auto FeatureParser::isValidFeature(const std::optional<std::vector<std::string>> &tokens)
    -> bool {
    const std::array<char, 3> allowedStrand = {'+', '-', '.'};

    return tokens.has_value() && tokens.value().size() == exptectedGffFileTokenCount &&
           std::ranges::find(allowedStrand, tokens.value()[strandTokenColumn][0]) !=
               allowedStrand.end();
}

}  // namespace annotation
