#include "SamFileUtility.hpp"

// Standard
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <iterator>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

// htslib
#include <htslib/hts.h>
#include <htslib/sam.h>

// seqan3
#include <seqan3/io/record.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/output.hpp>

// Internal
#include "Logger.hpp"

namespace SamFileUtility {

extern "C" {
using SamOrder = enum {
    Coordinate,
    QueryName,
    TagCoordinate,
    TagQueryName,
    MinHash,
    TemplateCoordinate
};
auto bam_sort_core_ext(SamOrder sam_order, char* sort_tag, int minimiser_kmer, bool try_rev,
                       bool no_squash, const char* fn, const char* prefix, const char* fnout,
                       const char* modeout, size_t _max_mem, int n_threads,
                       const htsFormat* in_fmt, const htsFormat* out_fmt, char* arg_list,
                       int no_pg, int write_index) -> int;
}

namespace {

using SamFilePtr = std::unique_ptr<samFile, decltype(&hts_close)>;
using SamHeaderPtr = std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)>;
using BamRecordPtr = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>;

[[nodiscard]] auto openSamFile(const fs::path& path) -> SamFilePtr {
    return SamFilePtr{sam_open(path.string().c_str(), "r"), hts_close};
}

auto fillReferenceData(const bam_hdr_t& header, SamFileInspection& inspection) -> void {
    inspection.referenceIDs.clear();
    inspection.referenceLengths.clear();

    if (header.n_targets <= 0) {
        return;
    }

    inspection.referenceIDs.resize(static_cast<std::size_t>(header.n_targets));
    inspection.referenceLengths.resize(static_cast<std::size_t>(header.n_targets));

    for (int32_t index = 0; index < header.n_targets; ++index) {
        const auto vectorIndex = static_cast<std::size_t>(index);
        inspection.referenceIDs[vectorIndex] = header.target_name[index];
        inspection.referenceLengths[vectorIndex] = header.target_len[index];
    }
}

auto uniqueSortToken() -> std::string {
    static std::atomic_uint64_t counter{0};
    const auto timestamp = std::chrono::steady_clock::now().time_since_epoch().count();
    return std::to_string(timestamp) + "." + std::to_string(counter++);
}

void removeSortTemporaryFiles(const fs::path& tempPrefix, const fs::path& tempOutput) {
    std::error_code ignoredError;
    if (!tempOutput.empty()) {
        fs::remove(tempOutput, ignoredError);
    }

    const auto parent = tempPrefix.parent_path().empty() ? fs::path{"."} : tempPrefix.parent_path();
    const auto prefixFilename = tempPrefix.filename().string();

    if (!fs::exists(parent)) {
        return;
    }

    for (const auto& entry : fs::directory_iterator(parent)) {
        if (!entry.is_regular_file()) {
            continue;
        }

        if (entry.path().filename().string().starts_with(prefixFilename)) {
            fs::remove(entry.path(), ignoredError);
        }
    }
}

}  // namespace

auto SamFileInspection::isReadable() const noexcept -> bool {
    switch (status) {
        case SamFileStatus::HeaderOnly:
        case SamFileStatus::HasRecords:
        case SamFileStatus::HeaderOnlyMissingEof:
        case SamFileStatus::HasRecordsMissingEof:
            return true;
        case SamFileStatus::Missing:
        case SamFileStatus::ZeroByte:
        case SamFileStatus::UnreadableOrTruncated:
            return false;
    }

    return false;
}

auto SamFileInspection::hasRecords() const noexcept -> bool {
    return status == SamFileStatus::HasRecords || status == SamFileStatus::HasRecordsMissingEof;
}

auto SamFileInspection::hasMissingEof() const noexcept -> bool {
    return status == SamFileStatus::HeaderOnlyMissingEof ||
           status == SamFileStatus::HasRecordsMissingEof;
}

auto SamFileInspection::hasReferenceHeader() const noexcept -> bool {
    return referenceIDs.size() == referenceLengths.size();
}

auto SamFileInspection::reference() const -> dataTypes::SamReference {
    return dataTypes::SamReference{referenceIDs, referenceLengths};
}

auto toString(const SamFileStatus status) noexcept -> std::string_view {
    switch (status) {
        case SamFileStatus::Missing:
            return "missing path";
        case SamFileStatus::ZeroByte:
            return "zero-byte file";
        case SamFileStatus::HeaderOnly:
            return "valid header-only BAM";
        case SamFileStatus::HasRecords:
            return "valid BAM with records";
        case SamFileStatus::HeaderOnlyMissingEof:
            return "readable header-only BAM missing BGZF EOF marker";
        case SamFileStatus::HasRecordsMissingEof:
            return "readable BAM with records missing BGZF EOF marker";
        case SamFileStatus::UnreadableOrTruncated:
            return "unreadable or truncated BAM";
    }

    return "unknown BAM status";
}

auto describe(const SamFileInspection& inspection) -> std::string {
    std::ostringstream out;
    out << "path=" << inspection.path << ", status=" << toString(inspection.status)
        << ", size=" << inspection.fileSize << ", records=" << inspection.recordCount
        << ", eof_check=" << inspection.eofCheck;

    if (!inspection.detail.empty()) {
        out << ", detail=" << inspection.detail;
    }

    return out.str();
}

auto inspect(const fs::path& path) -> SamFileInspection {
    SamFileInspection inspection{};
    inspection.path = path;

    if (!fs::exists(path)) {
        inspection.status = SamFileStatus::Missing;
        inspection.detail = "path does not exist";
        return inspection;
    }

    inspection.fileSize = fs::file_size(path);
    if (inspection.fileSize == 0) {
        inspection.status = SamFileStatus::ZeroByte;
        inspection.detail = "file is empty";
        return inspection;
    }

    {
        auto file = openSamFile(path);
        if (!file) {
            inspection.status = SamFileStatus::UnreadableOrTruncated;
            inspection.detail = "could not open file";
            return inspection;
        }

        inspection.eofCheck = hts_check_EOF(file.get());
    }

    auto file = openSamFile(path);
    if (!file) {
        inspection.status = SamFileStatus::UnreadableOrTruncated;
        inspection.detail = "could not reopen file after EOF check";
        return inspection;
    }

    SamHeaderPtr header{sam_hdr_read(file.get()), sam_hdr_destroy};
    if (!header) {
        inspection.status = SamFileStatus::UnreadableOrTruncated;
        inspection.detail = "could not read SAM/BAM header";
        return inspection;
    }

    fillReferenceData(*header, inspection);

    BamRecordPtr record{bam_init1(), bam_destroy1};
    if (!record) {
        inspection.status = SamFileStatus::UnreadableOrTruncated;
        inspection.detail = "could not allocate BAM record";
        return inspection;
    }

    int readResult = 0;
    while ((readResult = sam_read1(file.get(), header.get(), record.get())) >= 0) {
        ++inspection.recordCount;
    }

    if (readResult != -1) {
        inspection.status = SamFileStatus::UnreadableOrTruncated;
        inspection.detail = "read failed before clean EOF";
        return inspection;
    }

    if (inspection.eofCheck == 0) {
        inspection.status = inspection.recordCount == 0 ? SamFileStatus::HeaderOnlyMissingEof
                                                        : SamFileStatus::HasRecordsMissingEof;
        return inspection;
    }

    inspection.status =
        inspection.recordCount == 0 ? SamFileStatus::HeaderOnly : SamFileStatus::HasRecords;
    return inspection;
}

void writeHeaderOnlyFile(const fs::path& path, const dataTypes::SamReference& reference) {
    seqan3::sam_file_output out{path, reference.referenceIDs, reference.referenceLengths};
}

void sortByQueryName(const fs::path& inputPath, const fs::path& outputPath,
                     const size_t threadCount) {
    const auto inputInspection = inspect(inputPath);
    if (!inputInspection.isReadable()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
            "Could not sort unreadable alignments; ", describe(inputInspection),
            "; output=", outputPath);
    }

    if (inputInspection.hasMissingEof()) {
        Logger::log<LogLevel::WARNING>(
            "Alignment BAM is readable but missing the BGZF EOF marker; sorting will rewrite it: ",
            describe(inputInspection));
    }

    constexpr size_t SORT_DEFAULT_MEGS_PER_THREAD = 768;
    const size_t maxMem = SORT_DEFAULT_MEGS_PER_THREAD << 20;
    const htsFormat inFmt = {sequence_data, bam, {.major = 1, .minor = 6}, no_compression, 0, 0};
    const htsFormat outFmt = {sequence_data, bam, {.major = 1, .minor = 6}, no_compression, 0, 0};

    const auto token = uniqueSortToken();
    const fs::path tempPrefix =
        outputPath.parent_path() / (outputPath.filename().string() + ".sort_chunks." + token);
    const fs::path tempOutput =
        outputPath.parent_path() / (outputPath.filename().string() + ".sorting." + token + ".bam");

    char emptyStr[] = "";       // NOLINT
    const char wbStr[] = "wb";  // NOLINT

    const auto inputPathString = inputPath.string();
    const auto tempPrefixString = tempPrefix.string();
    const auto tempOutputString = tempOutput.string();
    const int sortThreadCount = static_cast<int>(std::max<size_t>(threadCount, 1));

    // NOLINTBEGIN
    int ret = bam_sort_core_ext(QueryName, emptyStr, 0, true, true, inputPathString.c_str(),
                                tempPrefixString.c_str(), tempOutputString.c_str(), wbStr, maxMem,
                                sortThreadCount, &inFmt, &outFmt, emptyStr, 1, 0);
    // NOLINTEND

    if (ret != 0) {
        removeSortTemporaryFiles(tempPrefix, tempOutput);
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
            "Could not sort alignments; input=", inputPath, "; output=", outputPath,
            "; temp_output=", tempOutput, "; input_status=", describe(inputInspection));
    }

    const auto sortedInspection = inspect(tempOutput);
    if (!sortedInspection.isReadable() || sortedInspection.hasMissingEof()) {
        removeSortTemporaryFiles(tempPrefix, tempOutput);
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
            "Sort produced invalid alignments; ", describe(sortedInspection),
            "; final_output=", outputPath);
    }

    try {
        fs::rename(tempOutput, outputPath);
    } catch (const std::exception& exception) {
        removeSortTemporaryFiles(tempPrefix, tempOutput);
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
            "Could not move sorted alignments into place; temp_output=", tempOutput,
            "; final_output=", outputPath, "; reason=", exception.what());
    }

    removeSortTemporaryFiles(tempPrefix, {});
}

auto countEntries(fs::path &path) -> std::size_t {
    seqan3::sam_file_input fin{path.string(), seqan3::fields<>{}};
    return std::ranges::distance(fin.begin(), fin.end());
}

auto hasAtLeastEntries(fs::path &path, size_t requiredCount) -> bool {
    seqan3::sam_file_input fin{path, seqan3::fields<>{}};

    size_t count{};

    for ([[maybe_unused]] auto const &entry : fin) {
        count++;

        if (count >= requiredCount) {
            return true;
        }
    }

    return false;
}

}  // namespace SamFileUtility
