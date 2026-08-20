#include "AnnotationHierarchyResolver.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <format>
#include <fstream>
#include <optional>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "AnnotationHierarchyError.hpp"
#include "GenomicFeature.hpp"
#include "GenomicStrand.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "TransparentStringMap.hpp"

namespace annotation {

namespace {

constexpr std::size_t kMaximumDiagnostics = 50;
constexpr std::size_t kMaximumRepairExamples = 10;

using FeatureIndex = std::size_t;
struct PendingEdge final {
    FeatureIndex childIndex{};
    std::string requestedParent;
    bool inferred{};
    std::optional<FeatureIndex> parentIndex;
};

struct AliasCandidate final {
    std::optional<FeatureIndex> first;
    std::optional<FeatureIndex> second;
};

using AliasCandidateMap =
    std::unordered_map<std::string, AliasCandidate, util::TransparentStringHash,
                       util::TransparentStringEqual>;

struct DiagnosticCollector final {
    std::size_t total{};
    std::vector<std::pair<std::size_t, std::string>> retained;

    auto add(std::size_t lineNumber, std::string message) -> void {
        ++total;
        if (retained.size() < kMaximumDiagnostics) {
            retained.emplace_back(lineNumber, std::move(message));
        }
    }
};

struct SourceMatch final {
    std::size_t lineNumber{};
    std::string featureType;
};

auto trim(std::string_view value) -> std::string_view {
    auto isSpace = [](const unsigned char character) { return std::isspace(character) != 0; };
    while (!value.empty() && isSpace(static_cast<unsigned char>(value.front()))) {
        value.remove_prefix(1);
    }
    while (!value.empty() && isSpace(static_cast<unsigned char>(value.back()))) {
        value.remove_suffix(1);
    }
    return value;
}

template <typename Callback>
auto forEachCommaSeparated(std::string_view values, Callback&& callback) -> void {
    while (!values.empty()) {
        const auto delimiter = values.find(',');
        const auto value = trim(values.substr(0, delimiter));
        if (!value.empty()) {
            callback(value);
        }
        if (delimiter == std::string_view::npos) {
            break;
        }
        values.remove_prefix(delimiter + 1);
    }
}

auto isChildFeatureType(std::string_view featureType) -> bool {
    return featureType == "exon" || featureType == "CDS" || featureType == "five_prime_UTR" ||
           featureType == "three_prime_UTR";
}

auto equivalent(const ParsedFeatureRecord& left, const ParsedFeatureRecord& right) -> bool {
    const auto& leftFeature = left.feature;
    const auto& rightFeature = right.feature;
    return leftFeature.getType() == rightFeature.getType() &&
           leftFeature.getGenomicRegion() == rightFeature.getGenomicRegion() &&
           leftFeature.getGeneName() == rightFeature.getGeneName() &&
           leftFeature.getAttributes() == rightFeature.getAttributes() &&
           left.explicitParentIds == right.explicitParentIds &&
           left.parentOrigin == right.parentOrigin;
}

auto structuredExonParent(std::string_view identifier) -> std::optional<std::string> {
    constexpr std::string_view prefix{"exon:"};
    if (!identifier.starts_with(prefix)) {
        return std::nullopt;
    }
    const auto parentStart = prefix.size();
    const auto delimiter = identifier.find(':', parentStart);
    if (delimiter == std::string_view::npos || delimiter == parentStart) {
        return std::nullopt;
    }
    return std::string{identifier.substr(parentStart, delimiter - parentStart)};
}

auto inferredParent(const ParsedFeatureRecord& record, FileType fileType)
    -> std::optional<std::string> {
    if (fileType != FileType::GFF || !isChildFeatureType(record.feature.getType())) {
        return std::nullopt;
    }

    const auto& attributes = record.feature.getAttributes();
    if (const auto transcript = attributes.find(std::string_view{"transcript_id"});
        transcript != attributes.end()) {
        const auto value = trim(transcript->second);
        if (!value.empty() && value != record.feature.getID()) {
            return std::string{value};
        }
    }

    if (record.feature.getType() == "exon") {
        return structuredExonParent(record.feature.getID());
    }
    return std::nullopt;
}

template <typename Callback>
auto forEachDeclaredAlias(const ParsedFeatureRecord& record, Callback&& callback) -> void {
    if (isChildFeatureType(record.feature.getType())) {
        return;
    }

    static constexpr std::array<std::string_view, 4> keys{"Name", "transcript_id", "gene_id",
                                                          "aliases"};
    const auto& attributes = record.feature.getAttributes();
    for (const auto key : keys) {
        const auto iterator = attributes.find(key);
        if (iterator == attributes.end()) {
            continue;
        }
        forEachCommaSeparated(iterator->second, callback);
    }

    const auto dbxref = attributes.find(std::string_view{"Dbxref"});
    if (dbxref == attributes.end()) {
        return;
    }
    forEachCommaSeparated(dbxref->second, [&](std::string_view value) {
        callback(value);
        const auto separator = value.find(':');
        if (separator != std::string_view::npos && separator + 1 < value.size()) {
            callback(value.substr(separator + 1));
        }
    });
}

auto registerAliasCandidate(AliasCandidate& candidate, FeatureIndex featureIndex) -> void {
    if (!candidate.first) {
        candidate.first = featureIndex;
    } else if (*candidate.first != featureIndex && !candidate.second) {
        candidate.second = featureIndex;
    }
}

auto splitColumns(std::string_view line) -> std::optional<std::array<std::string_view, 9>> {
    std::array<std::string_view, 9> columns{};
    for (std::size_t index = 0; index < columns.size(); ++index) {
        const auto delimiter = line.find('\t');
        if (delimiter == std::string_view::npos) {
            if (index + 1 != columns.size()) {
                return std::nullopt;
            }
            columns[index] = line;
            return columns;
        }
        columns[index] = line.substr(0, delimiter);
        line.remove_prefix(delimiter + 1);
    }
    return std::nullopt;
}

template <typename Callback>
auto forEachRawAttribute(std::string_view attributes, Callback&& callback) -> void {
    while (!attributes.empty()) {
        const auto delimiter = attributes.find(';');
        const auto field = trim(attributes.substr(0, delimiter));
        if (delimiter == std::string_view::npos) {
            attributes = {};
        } else {
            attributes.remove_prefix(delimiter + 1);
        }
        const auto assignment = field.find('=');
        if (assignment == std::string_view::npos) {
            continue;
        }
        callback(trim(field.substr(0, assignment)), trim(field.substr(assignment + 1)));
    }
}

auto findSourceMatches(const std::filesystem::path& path,
                       const std::unordered_set<std::string>& missingKeys)
    -> std::unordered_map<std::string, SourceMatch> {
    std::unordered_map<std::string, SourceMatch> matches;
    if (missingKeys.empty()) {
        return matches;
    }

    std::ifstream input{path};
    std::string line;
    std::size_t lineNumber{};
    while (std::getline(input, line)) {
        ++lineNumber;
        if (line.empty() || line.front() == '#') {
            continue;
        }
        const auto columns = splitColumns(line);
        if (!columns) {
            continue;
        }
        if (isChildFeatureType((*columns)[2])) {
            continue;
        }

        auto consider = [&](std::string_view value) {
            if (missingKeys.contains(std::string{value}) && !matches.contains(std::string{value})) {
                matches.emplace(std::string{value},
                                SourceMatch{lineNumber, std::string{(*columns)[2]}});
            }
        };
        forEachRawAttribute((*columns)[8], [&](std::string_view key, std::string_view value) {
            if (key == "ID" || key == "Name" || key == "transcript_id" || key == "gene_id" ||
                key == "aliases") {
                forEachCommaSeparated(value, consider);
            } else if (key == "Dbxref") {
                forEachCommaSeparated(value, [&](std::string_view token) {
                    consider(token);
                    const auto separator = token.find(':');
                    if (separator != std::string_view::npos && separator + 1 < token.size()) {
                        consider(token.substr(separator + 1));
                    }
                });
            }
        });
    }
    return matches;
}

auto formatFailure(const std::filesystem::path& path, const DiagnosticCollector& diagnostics)
    -> std::string {
    std::ostringstream message;
    message << "Annotation hierarchy validation failed for " << path.string() << " with "
            << diagnostics.total << " error" << (diagnostics.total == 1 ? "" : "s") << ".";
    for (const auto& [lineNumber, diagnostic] : diagnostics.retained) {
        message << "\n  Line " << lineNumber << ": " << diagnostic;
    }
    if (diagnostics.total > diagnostics.retained.size()) {
        message << "\n  ... " << (diagnostics.total - diagnostics.retained.size())
                << " additional errors omitted.";
    }
    message << "\nNo masking or alignment was performed.";
    return message.str();
}

}  // namespace

auto AnnotationHierarchyResolver::resolve(std::vector<ParsedFeatureRecord>&& records,
                                          const std::filesystem::path& featureFilePath,
                                          FileType fileType,
                                          const std::unordered_set<std::string>& includedFeatures)
    -> HierarchyResolutionResult {
    HierarchyResolutionResult result;

    // GTF uses transcript_id as its format-defined grouping key and commonly
    // repeats gene_id across records. Preserve the established GTF behavior;
    // strict repair is specific to GFF3's explicit ID/Parent hierarchy.
    if (fileType == FileType::GTF) {
        result.features.reserve(records.size());
        for (auto& record : records) {
            std::optional<std::string> parentId;
            if (!record.explicitParentIds.empty()) {
                parentId = record.explicitParentIds.front();
            }
            record.feature.groupID = std::move(parentId);
            result.features.emplace_back(std::move(record.feature));
        }
        return result;
    }

    DiagnosticCollector diagnostics;

    std::vector<ParsedFeatureRecord> uniqueRecords;
    uniqueRecords.reserve(records.size());
    std::unordered_map<std::string, FeatureIndex> identifierToIndex;
    identifierToIndex.reserve(records.size());

    for (auto& record : records) {
        const std::string& identifier = record.feature.getID();
        const auto existing = identifierToIndex.find(identifier);
        if (existing == identifierToIndex.end()) {
            const auto index = uniqueRecords.size();
            identifierToIndex.emplace(identifier, index);
            uniqueRecords.emplace_back(std::move(record));
            continue;
        }
        if (equivalent(uniqueRecords[existing->second], record)) {
            ++result.stats.equivalentDuplicates;
        } else {
            diagnostics.add(record.lineNumber,
                            std::format("feature ID '{}' duplicates line {} with conflicting data.",
                                        identifier, uniqueRecords[existing->second].lineNumber));
        }
    }

    std::vector<PendingEdge> edges;
    edges.reserve(uniqueRecords.size());
    AliasCandidateMap unresolvedAliases;

    for (FeatureIndex index = 0; index < uniqueRecords.size(); ++index) {
        auto& record = uniqueRecords[index];
        std::ranges::sort(record.explicitParentIds);
        record.explicitParentIds.erase(std::ranges::unique(record.explicitParentIds).begin(),
                                       record.explicitParentIds.end());
        if (record.explicitParentIds.size() > 1) {
            diagnostics.add(record.lineNumber,
                            std::format("feature '{}' declares multiple Parents; RNAnue requires a "
                                        "tree and cannot select one without information loss.",
                                        record.feature.getID()));
            continue;
        }

        std::optional<std::string> requestedParent;
        bool inferred = false;
        if (!record.explicitParentIds.empty()) {
            requestedParent = record.explicitParentIds.front();
        } else {
            requestedParent = inferredParent(record, fileType);
            inferred = requestedParent.has_value();
        }
        if (!requestedParent) {
            continue;
        }

        PendingEdge edge{.childIndex = index,
                         .requestedParent = std::move(*requestedParent),
                         .inferred = inferred,
                         .parentIndex = std::nullopt};
        if (const auto exact = identifierToIndex.find(edge.requestedParent);
            exact != identifierToIndex.end()) {
            edge.parentIndex = exact->second;
        } else {
            unresolvedAliases.try_emplace(edge.requestedParent);
        }
        edges.emplace_back(std::move(edge));
    }

    if (!unresolvedAliases.empty()) {
        for (FeatureIndex index = 0; index < uniqueRecords.size(); ++index) {
            forEachDeclaredAlias(uniqueRecords[index], [&](std::string_view alias) {
                const auto unresolved = unresolvedAliases.find(alias);
                if (unresolved != unresolvedAliases.end()) {
                    registerAliasCandidate(unresolved->second, index);
                }
            });
        }
    }

    std::unordered_set<std::string> missingKeys;
    for (auto& edge : edges) {
        if (edge.parentIndex) {
            continue;
        }
        const auto& candidate = unresolvedAliases.at(edge.requestedParent);
        if (candidate.first && !candidate.second) {
            edge.parentIndex = candidate.first;
            ++result.stats.aliasRepairs;
            if (result.stats.repairExamples.size() < kMaximumRepairExamples) {
                result.stats.repairExamples.emplace_back(std::format(
                    "feature '{}' Parent '{}' resolved to '{}' through a "
                    "declared alias",
                    uniqueRecords[edge.childIndex].feature.getID(), edge.requestedParent,
                    uniqueRecords[*candidate.first].feature.getID()));
            }
        } else if (candidate.second) {
            diagnostics.add(
                uniqueRecords[edge.childIndex].lineNumber,
                std::format("feature '{}' Parent '{}' is ambiguous between '{}' "
                            "(line {}) and '{}' "
                            "(line {}).",
                            uniqueRecords[edge.childIndex].feature.getID(), edge.requestedParent,
                            uniqueRecords[*candidate.first].feature.getID(),
                            uniqueRecords[*candidate.first].lineNumber,
                            uniqueRecords[*candidate.second].feature.getID(),
                            uniqueRecords[*candidate.second].lineNumber));
        } else {
            missingKeys.insert(edge.requestedParent);
        }
    }

    const auto sourceMatches = findSourceMatches(featureFilePath, missingKeys);
    for (const auto& edge : edges) {
        if (edge.parentIndex || !missingKeys.contains(edge.requestedParent)) {
            continue;
        }
        const auto& child = uniqueRecords[edge.childIndex];
        std::string guidance;
        if (const auto match = sourceMatches.find(edge.requestedParent);
            match != sourceMatches.end()) {
            const bool excluded =
                !includedFeatures.empty() && !includedFeatures.contains(match->second.featureType);
            guidance = excluded ? std::format(
                                      " A matching '{}' record exists at line {}; add '{}' to "
                                      "--featuretypes.",
                                      match->second.featureType, match->second.lineNumber,
                                      match->second.featureType)
                                : std::format(
                                      " A matching record exists at line {} but was not "
                                      "retained; "
                                      "check its reference ID and formatting.",
                                      match->second.lineNumber);
        }
        diagnostics.add(
            child.lineNumber,
            std::format("feature '{}' has unresolved {}Parent '{}'.{}", child.feature.getID(),
                        edge.inferred ? "inferred " : "", edge.requestedParent, guidance));
    }

    std::vector<std::optional<FeatureIndex>> parentByIndex(uniqueRecords.size());
    for (const auto& edge : edges) {
        if (!edge.parentIndex) {
            continue;
        }
        const auto childIndex = edge.childIndex;
        const auto parentIndex = *edge.parentIndex;
        const auto& child = uniqueRecords[childIndex];
        const auto& parent = uniqueRecords[parentIndex];
        if (childIndex == parentIndex) {
            diagnostics.add(child.lineNumber,
                            std::format("feature '{}' is its own Parent.", child.feature.getID()));
            continue;
        }

        const auto childRegion = child.feature.getGenomicRegion();
        const auto parentRegion = parent.feature.getGenomicRegion();
        if (childRegion.getReferenceIDIndex() != parentRegion.getReferenceIDIndex()) {
            diagnostics.add(child.lineNumber,
                            std::format("feature '{}' and Parent '{}' use different references.",
                                        child.feature.getID(), parent.feature.getID()));
            continue;
        }
        if (childRegion.getStrand() != dataTypes::GenomicStrand::NONE &&
            parentRegion.getStrand() != dataTypes::GenomicStrand::NONE &&
            childRegion.getStrand() != parentRegion.getStrand()) {
            diagnostics.add(child.lineNumber,
                            std::format("feature '{}' and Parent '{}' have incompatible strands.",
                                        child.feature.getID(), parent.feature.getID()));
            continue;
        }
        if (childRegion.getStart() < parentRegion.getStart() ||
            childRegion.getEnd() > parentRegion.getEnd()) {
            diagnostics.add(child.lineNumber,
                            std::format("feature '{}' lies outside Parent '{}' coordinates.",
                                        child.feature.getID(), parent.feature.getID()));
            continue;
        }
        parentByIndex[childIndex] = parentIndex;
        if (edge.inferred) {
            ++result.stats.inferredParents;
            if (result.stats.repairExamples.size() < kMaximumRepairExamples) {
                result.stats.repairExamples.emplace_back(
                    std::format("feature '{}' inferred Parent '{}'", child.feature.getID(),
                                parent.feature.getID()));
            }
        }
    }

    std::vector<unsigned char> visitState(uniqueRecords.size(), 0);
    for (FeatureIndex start = 0; start < uniqueRecords.size(); ++start) {
        if (visitState[start] != 0) {
            continue;
        }
        std::vector<FeatureIndex> path;
        FeatureIndex current = start;
        while (visitState[current] == 0) {
            visitState[current] = 1;
            path.push_back(current);
            if (!parentByIndex[current]) {
                break;
            }
            current = *parentByIndex[current];
        }
        if (visitState[current] == 1 && parentByIndex[current]) {
            diagnostics.add(uniqueRecords[current].lineNumber,
                            std::format("cycle detected involving feature '{}' and Parent '{}'.",
                                        uniqueRecords[current].feature.getID(),
                                        uniqueRecords[*parentByIndex[current]].feature.getID()));
        }
        for (const auto index : path) {
            visitState[index] = 2;
        }
    }

    if (diagnostics.total != 0) {
        throw AnnotationHierarchyError{formatFailure(featureFilePath, diagnostics)};
    }

    // Materialize every canonical parent ID before moving any feature. A parent
    // can precede its children, so reading its ID during the move loop would
    // otherwise observe a moved-from string.
    for (FeatureIndex index = 0; index < uniqueRecords.size(); ++index) {
        std::optional<std::string> parentId;
        if (parentByIndex[index]) {
            parentId = uniqueRecords[*parentByIndex[index]].feature.getID();
        }
        uniqueRecords[index].feature.groupID = std::move(parentId);
    }

    result.features.reserve(uniqueRecords.size());
    for (auto& record : uniqueRecords) {
        result.features.emplace_back(std::move(record.feature));
    }

    const auto repairCount = result.stats.inferredParents + result.stats.aliasRepairs;
    if (repairCount != 0 || result.stats.equivalentDuplicates != 0) {
        Logger::log<LogLevel::WARNING>(
            std::format("Annotation hierarchy repaired: {} inferred Parents, {} "
                        "alias resolutions, {} "
                        "equivalent duplicates removed.",
                        result.stats.inferredParents, result.stats.aliasRepairs,
                        result.stats.equivalentDuplicates));
        for (const auto& example : result.stats.repairExamples) {
            Logger::log<LogLevel::DEBUG>(example);
        }
    }

    return result;
}

}  // namespace annotation
