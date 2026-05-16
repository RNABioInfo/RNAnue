#pragma once

// Standard
#include <cstddef>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

// Internal
#include "SamReference.hpp"

namespace SamFileUtility {

namespace fs = std::filesystem;

enum class SamFileStatus : std::uint8_t {
    Missing,
    ZeroByte,
    HeaderOnly,
    HasRecords,
    HeaderOnlyMissingEof,
    HasRecordsMissingEof,
    UnreadableOrTruncated
};

struct SamFileInspection {
    fs::path path;
    SamFileStatus status{SamFileStatus::Missing};
    std::uintmax_t fileSize{0};
    std::size_t recordCount{0};
    int eofCheck{-1};
    std::deque<std::string> referenceIDs;
    std::vector<size_t> referenceLengths;
    std::string detail;

    [[nodiscard]] auto isReadable() const noexcept -> bool;
    [[nodiscard]] auto hasRecords() const noexcept -> bool;
    [[nodiscard]] auto hasMissingEof() const noexcept -> bool;
    [[nodiscard]] auto hasReferenceHeader() const noexcept -> bool;
    [[nodiscard]] auto reference() const -> dataTypes::SamReference;
};

[[nodiscard]] auto toString(SamFileStatus status) noexcept -> std::string_view;
[[nodiscard]] auto describe(const SamFileInspection& inspection) -> std::string;

[[nodiscard]] auto inspect(const fs::path& path) -> SamFileInspection;

void writeHeaderOnlyFile(const fs::path& path, const dataTypes::SamReference& reference);

void sortByQueryName(const fs::path& inputPath, const fs::path& outputPath, size_t threadCount);

auto countEntries(fs::path& path) -> std::size_t;

auto hasAtLeastEntries(fs::path& path, size_t requiredCount) -> bool;

}  // namespace SamFileUtility
