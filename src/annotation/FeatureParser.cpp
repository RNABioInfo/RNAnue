// FeatureParser.cpp
#include "FeatureParser.hpp"

// Standard
#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

// Internal
#include "AnnotationHierarchyResolver.hpp"
#include "Constants.hpp"
#include "FeatureGrouper.hpp"
#include "FileType.hpp"
#include "GenomicFeature.hpp"
#include "GenomicFeatureGroup.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "ReferenceIndexMapping.hpp"

namespace annotation {

using namespace constants::annotation;

FeatureParser::FeatureParser(IncludedFeatureSet includedFeaturesParam,
                             std::optional<std::string> featureIdFlagParam)
    : includedFeatures(std::move(includedFeaturesParam)),
      featureIDFlag(std::move(featureIdFlagParam)) {}

auto FeatureParser::parseFlatMap(const fs::path& featureFilePath,
                                 const ReferenceIndexMapping& referenceToIndex) const
    -> FeatureMap {
    return parseFlatAndGrouped(featureFilePath, referenceToIndex).flatByChromosomeIndex;
}

auto FeatureParser::parseGroupedByParentID(const fs::path& featureFilePath,
                                           const ReferenceIndexMapping& referenceToIndex) const
    -> ParentIDToFeatureGroupMap {
    return parseFlatAndGrouped(featureFilePath, referenceToIndex).groupedByParentID;
}

auto FeatureParser::parseGroupedByHierarchy(const fs::path& featureFilePath,
                                            const ReferenceIndexMapping& referenceToIndex) const
    -> ParentIDToFeatureGroupMap {
    static_assert(columnCount == expectedAnnotationFileTokenCount);

    const auto fileType = getFileType(featureFilePath);

    ParseSettings settings{
        .fileType = fileType,
        .idKey = std::string_view{featureIDFlag.value_or(fileType.defaultIDKey())},
        .keys =
            AttributeKeys{
                .parentKey = std::string_view{fileType.defaultGroupKey()},
                .geneNameKey = std::string_view{FileType::defaultGeneNameKey()},
            },
    };

    ParentIDToFeatureGroupMap result{};
    const auto stats = scanFile(FileScanInput{
        .featureFilePath = featureFilePath,
        .settings = settings,
        .referenceToIndex = referenceToIndex,
        .flatMap = nullptr,
        .groupMap = &result,
        .groupingMode = GroupingMode::Hierarchy,
    });

    const auto includedText = buildIncludedFeatureTypesText(includedFeatures);
    if (stats.parsedCount == 0) {
        Logger::log<LogLevel::WARNING>(
            std::format("No features parsed from file {}. Check your feature file and included "
                        "features flag: {}",
                        featureFilePath.string(), includedText));
    } else {
        const auto groupInfo =
            stats.parentIDs.empty()
                ? std::string{}
                : (std::string{" Found "} + std::to_string(stats.parentIDs.size()) +
                   " parent groups.");

        Logger::log(
            std::format("Parsed {} features of type: {}.{}", stats.parsedCount, includedText,
                        groupInfo));
    }

    return result;
}

auto FeatureParser::parseFlatAndGrouped(const fs::path& featureFilePath,
                                        const ReferenceIndexMapping& referenceToIndex) const
    -> Results {
    static_assert(columnCount == expectedAnnotationFileTokenCount);

    const auto fileType = getFileType(featureFilePath);

    ParseSettings settings{
        .fileType = fileType,
        .idKey = std::string_view{featureIDFlag.value_or(fileType.defaultIDKey())},
        .keys =
            AttributeKeys{
                .parentKey = std::string_view{fileType.defaultGroupKey()},
                .geneNameKey = std::string_view{FileType::defaultGeneNameKey()},
            },
    };

    Results result{};
    const auto stats = scanFile(FileScanInput{
        .featureFilePath = featureFilePath,
        .settings = settings,
        .referenceToIndex = referenceToIndex,
        .flatMap = &result.flatByChromosomeIndex,
        .groupMap = &result.groupedByParentID,
        .groupingMode = GroupingMode::DirectParentID,
    });

    const auto includedText = buildIncludedFeatureTypesText(includedFeatures);

    if (stats.parsedCount == 0) {
        Logger::log<LogLevel::WARNING>(
            std::format("No features parsed from file {}. Check your feature file and included "
                        "features flag: {}",
                        featureFilePath.string(), includedText));
        return result;
    }

    const auto groupInfo =
        stats.parentIDs.empty()
            ? std::string{}
            : (std::string{" Found "} + std::to_string(stats.parentIDs.size()) + " parent groups.");

    const auto msg = std::format("Parsed {} features of type: {}.{}", stats.parsedCount,
                                 includedText, groupInfo);
    Logger::log(msg);

    return result;
}

auto FeatureParser::getFileType(const fs::path& featureFilePath) -> FileType {
    std::ifstream fileStream(featureFilePath);
    if (!fileStream) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not open annotation file: {}",
                                                            featureFilePath.string());
    }

    std::string lineData;
    while (std::getline(fileStream, lineData)) {
        if (!lineData.empty()) {
            break;
        }
    }

    if (lineData.starts_with("##gff-version")) {
        return FileType{FileType::GFF};
    }
    if (lineData.starts_with("##gtf-version")) {
        return FileType{FileType::GTF};
    }

    Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
        "Annotation file type not supported. First non-empty line of {} does not contain "
        "##gff-version or ##gtf-version",
        featureFilePath.string());

    std::unreachable();
}

auto FeatureParser::shouldSkipLine(std::string_view line) -> bool {
    return line.empty() || line.front() == '#';
}

auto FeatureParser::trim(std::string_view value) -> std::string_view {
    auto isSpace = [](unsigned char chr) { return std::isspace(chr) != 0; };

    while (!value.empty() && isSpace(static_cast<unsigned char>(value.front()))) {
        value.remove_prefix(1);
    }
    while (!value.empty() && isSpace(static_cast<unsigned char>(value.back()))) {
        value.remove_suffix(1);
    }
    return value;
}

auto FeatureParser::tryParseInt(std::string_view value) -> std::optional<int> {
    int parsed{};
    const char* begPtr = value.data();
    const char* endPtr = value.data() + value.size();

    const auto res = std::from_chars(begPtr, endPtr, parsed);
    if (res.ec != std::errc{} || res.ptr != endPtr) {
        return std::nullopt;
    }
    return parsed;
}

auto FeatureParser::tryParseCoordinates(CoordinateTokens tokens)
    -> std::optional<std::pair<int, int>> {
    const auto startOne = tryParseInt(tokens.startToken);
    const auto endOne = tryParseInt(tokens.endToken);

    if (!startOne || !endOne || *startOne <= 0 || *endOne <= 0) {
        return std::nullopt;
    }

    return std::pair<int, int>{*startOne - 1, *endOne};
}

auto FeatureParser::tryParseStrand(std::string_view token) -> std::optional<char> {
    if (token.empty()) {
        return std::nullopt;
    }

    const char chr = token.front();
    if (chr == '+' || chr == '-' || chr == '.') {
        return chr;
    }

    return std::nullopt;
}

auto FeatureParser::trySplitColumns(std::string_view line, const char delimiter)
    -> std::optional<ColumnViews> {
    ColumnViews out{};

    for (std::size_t idxCol = 0; idxCol < columnCount; ++idxCol) {
        const auto posDelim = line.find(delimiter);

        if (posDelim == std::string_view::npos) {
            if (idxCol + 1 != columnCount) {
                return std::nullopt;
            }
            out.cols[idxCol] = line;
            break;
        }

        out.cols[idxCol] = line.substr(0, posDelim);
        line.remove_prefix(posDelim + 1);
    }

    if (out.cols.back().find(delimiter) != std::string_view::npos) {
        return std::nullopt;
    }

    return out;
}

auto FeatureParser::normalizeAttributeValue(const FileType fileType, std::string_view value)
    -> std::string {
    value = trim(value);

    if (fileType != FileType::GTF) {
        return std::string{value};
    }

    std::string out;
    out.reserve(value.size());

    for (char chr : value) {
        if (chr != '"') {
            out.push_back(chr);
        }
    }

    while (!out.empty()) {
        const auto lastChr = static_cast<unsigned char>(out.back());
        if (out.back() == ';' || std::isspace(lastChr) != 0) {
            out.pop_back();
            continue;
        }
        break;
    }

    const auto trimmed = trim(std::string_view{out});
    if (trimmed.size() == out.size()) {
        return out;
    }
    return std::string{trimmed};
}

auto FeatureParser::extractAllAttributes(AttributeParseInput const& input) -> AttributeMap {
    const char fieldDelim = FileType::attributeDelimiter();
    const char assignChar = input.fileType.attributeAssignment();

    AttributeMap attributes;
    const auto estimatedCount = estimateAttributeCount(input.attributes, fieldDelim);
    attributes.reserve(std::max(kMinAttributeReserve, estimatedCount));

    std::string_view remaining = input.attributes;

    while (!remaining.empty()) {
        const auto delimPos = remaining.find(fieldDelim);
        const auto rawField =
            remaining.substr(0, delimPos == std::string_view::npos ? remaining.size() : delimPos);
        const auto fieldView = trim(rawField);

        if (delimPos == std::string_view::npos) {
            remaining = {};
        } else {
            remaining.remove_prefix(delimPos + 1);
        }

        if (fieldView.empty()) {
            continue;
        }

        const auto assignPos = fieldView.find(assignChar);
        if (assignPos == std::string_view::npos) {
            continue;
        }

        const auto keyView = trim(fieldView.substr(0, assignPos));
        const auto valView = trim(fieldView.substr(assignPos + 1));
        if (keyView.empty() || valView.empty()) {
            continue;
        }

        attributes.try_emplace(std::string{keyView}, std::string{valView});
    }

    return attributes;
}

auto FeatureParser::tryParseFeature(const ColumnViews& columns, const ParseSettings& settings,
                                    const ReferenceIndexMapping& referenceToIndex,
                                    LookupBuffers& buffers, std::size_t lineNumber) const
    -> std::optional<ParsedFeature> {
    const auto& cols = columns.cols;

    const std::string_view refView = cols[0];
    const std::string_view typView = cols[2];

    if (!includedFeatures.empty()) {
        buffers.featureType.assign(typView.data(), typView.size());
        if (!includedFeatures.contains(buffers.featureType)) {
            return std::nullopt;
        }
    }

    buffers.referenceID.assign(refView.data(), refView.size());
    const auto refIdxOpt = referenceToIndex.findIndex(buffers.referenceID);
    if (!refIdxOpt) {
        Logger::log<LogLevel::WARNING>(
            std::format("Feature reference id not found: {}. Ensure the annotation and reference "
                        "genome use the same reference IDs.",
                        buffers.referenceID));

        return std::nullopt;
    }

    const int refIdx = *refIdxOpt;

    const auto coordOpt =
        tryParseCoordinates(CoordinateTokens{.startToken = cols[3], .endToken = cols[4]});
    if (!coordOpt) {
        Logger::log<LogLevel::WARNING>("Could not parse coordinates for feature.");
        return std::nullopt;
    }

    const auto strandOpt = tryParseStrand(cols[strandTokenColumn]);
    if (!strandOpt) {
        Logger::log<LogLevel::WARNING>("Could not parse strand for feature.");
        return std::nullopt;
    }

    const auto attrs = extractAttributes(AttributeParseInput{
        .fileType = settings.fileType,
        .attributes = cols[8],
        .idKey = settings.idKey,
        .keys = settings.keys,
    });

    if (!attrs.identifier) {
        Logger::log<LogLevel::WARNING>(std::format("Could not parse identifier for feature."));
        return std::nullopt;
    }

    GenomicRegion region{
        refIdx,
        {.startPosition = coordOpt->first, .endPosition = coordOpt->second},
        getGenomicStrand(*strandOpt),
    };

    return ParsedFeature{
        .feature = GenomicFeature{std::string(typView), region, *attrs.identifier, std::nullopt,
                                  attrs.geneName, attrs.attributes},
        .explicitParentIds = attrs.parentIDs,
        .lineNumber = lineNumber,
        .parentOrigin = attrs.parentIDs.empty()
                            ? ParsedParentOrigin::None
                            : (settings.fileType == FileType::GTF
                                   ? ParsedParentOrigin::GtfTranscriptId
                                   : ParsedParentOrigin::GffParent),
    };
}

auto FeatureParser::scanFile(FileScanInput input) const -> ScanStats {
    std::ifstream fileStream(input.featureFilePath);
    if (!fileStream) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not open annotation file: {}",
                                                            input.featureFilePath.string());
    }

    ScanStats scanStats{};
    LookupBuffers buffers{};

    std::vector<ParsedFeatureRecord> parsedRecords;
    constexpr std::size_t initialReserve = 1'024;
    parsedRecords.reserve(initialReserve);

    std::string lineData;
    std::size_t lineIndex = 0;

    while (std::getline(fileStream, lineData)) {
        ++lineIndex;

        const std::string_view lineView{lineData};
        if (shouldSkipLine(lineView)) {
            continue;
        }

        const auto columnsOpt = trySplitColumns(lineView, '\t');
        if (!columnsOpt) {
            Logger::log<LogLevel::WARNING>(
                "Skipping malformed feature at line {} in {}: expected {} columns", lineIndex,
                input.featureFilePath.string(), expectedAnnotationFileTokenCount);
            continue;
        }

        auto parsedOpt = tryParseFeature(*columnsOpt, input.settings, input.referenceToIndex,
                                         buffers, lineIndex);
        if (!parsedOpt) {
            continue;
        }

        parsedRecords.emplace_back(std::move(*parsedOpt));
    }

    auto resolution = AnnotationHierarchyResolver::resolve(
        std::move(parsedRecords), input.featureFilePath, input.settings.fileType, includedFeatures);
    scanStats.parsedCount = resolution.features.size();

    if (input.flatMap != nullptr) {
        for (const auto& feature : resolution.features) {
            auto& vecRef = (*input.flatMap)[feature.getGenomicRegion().getReferenceIDIndex()];
            vecRef.emplace_back(feature);
        }
    }

    if (input.groupMap != nullptr) {
        auto groupingResult =
            input.groupingMode == GroupingMode::Hierarchy
                ? (input.settings.fileType == FileType::GFF
                       ? annotation::FeatureGrouper::groupByValidatedHierarchy(
                             std::move(resolution.features))
                       : annotation::FeatureGrouper::groupByHierarchy(
                             std::move(resolution.features)))
                : annotation::FeatureGrouper::groupByDirectParentID(std::move(resolution.features));

        *input.groupMap = std::move(groupingResult.groups);
        scanStats.parentIDs = std::move(groupingResult.groupKeys);
    }

    return scanStats;
}

auto FeatureParser::buildIncludedFeatureTypesText(const IncludedFeatureSet& featureSet)
    -> std::string {
    std::string text;
    constexpr int capacity = 64;
    text.reserve(capacity);

    for (const auto& item : featureSet) {
        if (!text.empty()) {
            text.append(", ");
        }
        text.append(item);
    }

    return text;
}
auto FeatureParser::estimateAttributeCount(std::string_view attributes, char fieldDelim)
    -> std::size_t {
    if (attributes.empty()) {
        return std::size_t{0};
    }

    std::size_t delimCount{0};
    for (const char chr : attributes) {
        delimCount += static_cast<std::size_t>(chr == fieldDelim);
    }
    return delimCount + std::size_t{1};
}
auto FeatureParser::popRaw(AttributeMap& attributes, std::string_view keyView)
    -> std::optional<std::string> {
    const auto iterator = attributes.find(keyView);
    if (iterator == attributes.end()) {
        return std::nullopt;
    }

    auto value = std::move(iterator->second);
    attributes.erase(iterator);
    return value;
}
auto FeatureParser::popNormalized(AttributeMap& attributes, FileType fileType,
                                  std::string_view keyView) -> std::optional<std::string> {
    auto rawValue = popRaw(attributes, keyView);
    if (!rawValue) {
        return std::nullopt;
    }
    return normalizeAttributeValue(fileType, std::string_view{*rawValue});
}
auto FeatureParser::popParentIds(AttributeMap& attributes, FileType fileType,
                                 std::string_view parentKeyView) -> std::vector<std::string> {
    auto normalized = popNormalized(attributes, fileType, parentKeyView);
    if (!normalized) {
        return {};
    }

    std::vector<std::string> parentIds;
    std::string_view remaining{*normalized};
    while (!remaining.empty()) {
        const auto commaPosition = remaining.find(',');
        const auto value = trim(remaining.substr(0, commaPosition));
        if (!value.empty()) {
            parentIds.emplace_back(value);
        }
        if (commaPosition == std::string_view::npos) {
            break;
        }
        remaining.remove_prefix(commaPosition + 1);
    }
    return parentIds;
}

auto FeatureParser::popSpecificAttributes(AttributeMap& attributes,
                                          AttributeParseInput const& input) -> ExtractedAttributes {
    ExtractedAttributes extracted;
    extracted.identifier = popNormalized(attributes, input.fileType, input.idKey);
    extracted.parentIDs = popParentIds(attributes, input.fileType, input.keys.parentKey);
    extracted.geneName = popNormalized(attributes, input.fileType, input.keys.geneNameKey);
    extracted.attributes = std::move(attributes);
    return extracted;
}
auto FeatureParser::extractAttributes(AttributeParseInput input) -> ExtractedAttributes {
    auto attributes = extractAllAttributes(input);
    return popSpecificAttributes(attributes, input);
}
}  // namespace annotation
