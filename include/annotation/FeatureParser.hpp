// FeatureParser.hpp
#pragma once

// Standard
#include <array>
#include <cctype>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>

// Internal
#include "Constants.hpp"
#include "FileType.hpp"
#include "GenomicFeature.hpp"
#include "GenomicFeatureGroup.hpp"
#include "ReferenceIndexMapping.hpp"
#include "TransparentStringMap.hpp"

namespace annotation {

namespace fs = std::filesystem;

using namespace dataTypes;

using IncludedFeatureSet = std::unordered_set<std::string>;

class FeatureParser {
   public:
    explicit FeatureParser(IncludedFeatureSet includedFeatures,
                           std::optional<std::string> featureIdFlag);

    FeatureParser(const FeatureParser&) = default;
    FeatureParser(FeatureParser&&) = delete;
    auto operator=(const FeatureParser&) -> FeatureParser& = delete;
    auto operator=(FeatureParser&&) -> FeatureParser& = delete;
    ~FeatureParser() = default;

    [[nodiscard]] static auto defaultAll() noexcept { return FeatureParser{{}, std::nullopt}; }
    [[nodiscard]] static auto defaultAllTranscript() -> FeatureParser {
        return FeatureParser{constants::annotation::defaultAllTranscriptTypes, std::nullopt};
    }

    using ResultFlat = dataTypes::FeatureMap;
    using ResultGrouped = dataTypes::ParentIDToFeatureGroupMap;

    struct Results {
        ResultFlat flatByChromosomeIndex;
        ResultGrouped groupedByParentID;
    };

    [[nodiscard]] auto parseFlatMap(const fs::path& featureFilePath,
                                    const ReferenceIndexMapping& referenceToIndex) const
        -> ResultFlat;

    [[nodiscard]] auto parseGroupedByParentID(const fs::path& featureFilePath,
                                              const ReferenceIndexMapping& referenceToIndex) const
        -> ResultGrouped;

    [[nodiscard]] auto parseFlatAndGrouped(const fs::path& featureFilePath,
                                           const ReferenceIndexMapping& referenceToIndex) const
        -> Results;

   private:
    IncludedFeatureSet includedFeatures;
    std::optional<std::string> featureIDFlag;

    static constexpr std::size_t columnCount = 9;

    struct AttributeKeys {
        std::string_view parentKey;
        std::string_view geneNameKey;
    };

    struct ParseSettings {
        FileType fileType;
        std::string_view idKey;
        AttributeKeys keys;
    };

    struct ColumnViews {
        std::array<std::string_view, columnCount> cols;
    };

    struct CoordinateTokens {
        std::string_view startToken;
        std::string_view endToken;
    };

    struct AttributeParseInput {
        FileType fileType;
        std::string_view attributes;
        std::string_view idKey;
        AttributeKeys keys;
    };

    struct ExtractedAttributes {
        std::optional<std::string> identifier;
        std::optional<std::string> parentID;
        std::optional<std::string> geneName;

        using AttributeMap = util::TransparentStringMap;
        AttributeMap attributes;
    };

    using ParsedFeature = GenomicFeature;

    struct LookupBuffers {
        std::string referenceID;
        std::string featureType;
    };

    struct ScanStats {
        std::size_t parsedCount{};
        std::unordered_set<std::string> parentIDs;
        std::unordered_set<std::string> nonAdjacentParentsWarned;
    };

    struct FileScanInput {
        const fs::path& featureFilePath;
        ParseSettings settings;
        const ReferenceIndexMapping& referenceToIndex;
        dataTypes::FeatureMap* flatMap;
        dataTypes::ParentIDToFeatureGroupMap* groupMap;
    };

    using AttributeMap = ExtractedAttributes::AttributeMap;

    static constexpr auto kMinAttributeReserve = std::size_t{8};

    [[nodiscard]] static auto getFileType(const fs::path& featureFilePath) -> FileType;

    static auto shouldSkipLine(std::string_view line) -> bool;
    static auto trim(std::string_view value) -> std::string_view;

    static auto tryParseInt(std::string_view value) -> std::optional<int>;
    static auto tryParseCoordinates(CoordinateTokens tokens) -> std::optional<std::pair<int, int>>;
    static auto tryParseStrand(std::string_view token) -> std::optional<char>;

    static auto trySplitColumns(std::string_view line, char delimiter)
        -> std::optional<ColumnViews>;

    static auto normalizeAttributeValue(FileType fileType, std::string_view value) -> std::string;
    static auto estimateAttributeCount(std::string_view attributes, char fieldDelim) -> std::size_t;
    static auto extractAllAttributes(AttributeParseInput const& input) -> AttributeMap;
    static auto popRaw(AttributeMap& attributes, std::string_view keyView)
        -> std::optional<std::string>;

    static auto popNormalized(AttributeMap& attributes, FileType fileType, std::string_view keyView)
        -> std::optional<std::string>;
    static auto popParentId(AttributeMap& attributes, FileType fileType,
                            std::string_view parentKeyView) -> std::optional<std::string>;
    static auto popSpecificAttributes(AttributeMap& attributes, AttributeParseInput const& input)
        -> ExtractedAttributes;
    static auto extractAttributes(AttributeParseInput input) -> ExtractedAttributes;

    [[nodiscard]] auto tryParseFeature(const ColumnViews& columns, const ParseSettings& settings,
                                       const ReferenceIndexMapping& referenceToIndex,
                                       LookupBuffers& buffers) const
        -> std::optional<ParsedFeature>;

    [[nodiscard]] auto scanFile(FileScanInput input) const -> ScanStats;

    static auto buildIncludedFeatureTypesText(const IncludedFeatureSet& featureSet) -> std::string;
};

}  // namespace annotation
