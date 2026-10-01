#include "SegemehlAligner.hpp"

// Standard
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

// segmehl
#include "AnnotationFilePicker.hpp"
#include "ReferenceGenomeFilePicker.hpp"
extern "C" {
#include "segemehl.h"
}

// Internal
#include "Logger.hpp"
#include "ReferenceGenome.hpp"
#include "SamFileUtility.hpp"
#include "SequenceFileUtility.hpp"

namespace pipelines::align {

auto SegemehlAligner::findIndex(const fs::path &referenceGenomePath) const
    -> std::optional<fs::path> {
    // Check if index exists in the same directory as the reference genome
    fs::path indexFileName = referenceGenomePath.filename().replace_extension(".idx");
    fs::path indexFilePath = referenceGenomePath.parent_path() / indexFileName;

    Logger::log("Searching for reference index at: ", indexFilePath);

    if (fs::exists(indexFilePath)) {
        return indexFilePath;
    }

    // Check if index exists in the output directory
    indexFilePath = parameters.outputDir / indexFileName;

    Logger::log("Searching for reference index at: ", indexFilePath);

    if (fs::exists(indexFilePath)) {
        return indexFilePath;
    }

    Logger::log("Did not find reference index");

    return std::nullopt;
}

void SegemehlAligner::buildIndex() {
    if (parameters.alignmentIndex) {
        indexPath = fs::absolute(*parameters.alignmentIndex);
        Logger::log("Using supplied segemehl index: ", indexPath);
        return;
    }
    fs::path referencePath = parameters.referenceGenome;
    size_t const threads = parameters.threadCount;

    const auto indexFilePath = findIndex(referencePath);

    if (indexFilePath.has_value()) {
        Logger::log("Existing index found: ", indexFilePath);
        indexPath = *indexFilePath;
        return;
    }

    // Index file is written to same location as reference genome
    indexPath = referencePath.parent_path() / referencePath.filename().replace_extension(".idx");

    Logger::log("Building index");
    std::vector<std::string> args = {"-x", indexPath.string(),     "-d", referencePath.string(),
                                     "-t", std::to_string(threads)};

    auto c_args = convertToCStrings(args);

    int result = segemehl(static_cast<int>(c_args.size()) - 1, c_args.data());

    if (result != 0) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not create index for: ",
                                                            referencePath);
    }
}

[[nodiscard]] auto SegemehlAligner::threadsAdaptedToEntries(const fs::path &inputPath) const
    -> size_t {
    const bool hasSufficientEntries =
        SequenceFileUtility::hasAtLeastEntries(inputPath, parameters.threadCount);

    size_t threads = parameters.threadCount;

    if (not hasSufficientEntries) {
        threads = SequenceFileUtility::countEntries(inputPath) > 0 ? 1 : 0;
    }

    return threads;
}

[[nodiscard]] auto SegemehlAligner::getGeneralAlignmentArgs(size_t threadCount) const
    -> std::vector<std::string> {
    return {"-b", "-S",
            "-A", std::to_string(parameters.accuracy),
            "-U", std::to_string(parameters.minimumFragmentScore),
            "-W", std::to_string(parameters.minimumSpliceCoverage),
            "-Z", std::to_string(parameters.minimumFragmentLength),
            "-t", std::to_string(threadCount),
            "-m", std::to_string(parameters.minLengthThreshold),
            "-i", indexPath.string(),
            "-d", parameters.referenceGenome.string()};
}

auto SegemehlAligner::referenceFromGenome() const -> dataTypes::SamReference {
    const ReferenceGenome referenceGenome{parameters.referenceGenome};

    auto referenceIDs = referenceGenome.getReferenceIndexMapping().sortedReferenceIDs();
    std::vector<size_t> referenceLengths;
    referenceLengths.reserve(referenceIDs.size());

    for (size_t index = 0; index < referenceIDs.size(); ++index) {
        referenceLengths.push_back(referenceGenome.getSequence(static_cast<int>(index)).size());
    }

    return dataTypes::SamReference{std::move(referenceIDs), std::move(referenceLengths)};
}

void SegemehlAligner::writeEmptyAlignments(const fs::path &alignmentsOutPath,
                                           const fs::path &emptyInputPath) const {
    Logger::log("File has no entries: ", emptyInputPath,
                "; writing header-only alignments: ", alignmentsOutPath);
    SamFileUtility::writeHeaderOnlyFile(alignmentsOutPath, referenceFromGenome());
}

void SegemehlAligner::runSegemehlAlignment(std::vector<std::string> args,
                                           const fs::path &outputPath,
                                           const std::string &errorMessage) {
    constexpr size_t VALIDATION_ATTEMPTS = 6;
    constexpr size_t INITIAL_RETRY_DELAY_MS = 500;
    constexpr size_t ALIGNMENT_ATTEMPTS = 2;

    for (size_t attempt = 1; attempt <= ALIGNMENT_ATTEMPTS; ++attempt) {
        auto c_args = convertToCStrings(args);

        const int result = segemehl(static_cast<int>(c_args.size()) - 1, c_args.data());

        if (result != 0) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(errorMessage);
        }

        const auto inspection = SamFileUtility::inspectWithRetries(outputPath, VALIDATION_ATTEMPTS,
                                                                   INITIAL_RETRY_DELAY_MS);
        if (inspection.isReadable()) {
            if (inspection.hasMissingEof()) {
                Logger::log<LogLevel::WARNING>(
                    "Alignment output is readable but missing the BGZF EOF marker; sorting will "
                    "repair it: ",
                    SamFileUtility::describe(inspection));
            }
            return;
        }

        if (attempt == ALIGNMENT_ATTEMPTS) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
                errorMessage,
                "; output failed validation after retry: ", SamFileUtility::describe(inspection));
        }

        Logger::log<LogLevel::WARNING>(
            "Alignment output failed validation after segemehl; regenerating once: ",
            SamFileUtility::describe(inspection));

        std::error_code ignoredError;
        fs::remove(outputPath, ignoredError);
    }
}

void SegemehlAligner::alignSingleReads(const fs::path &queryFastqInPath,
                                       const fs::path &alignmentsFastqOutPath) const {
    const size_t threads = threadsAdaptedToEntries(queryFastqInPath);

    if (threads == 0) {
        writeEmptyAlignments(alignmentsFastqOutPath, queryFastqInPath);
        return;
    }

    auto args = getGeneralAlignmentArgs(threads);

    args.insert(args.end(),
                {"-q", queryFastqInPath.string(), "-o", alignmentsFastqOutPath.string(), "-H",
                 std::to_string(static_cast<int>(!parameters.multimapAlignments))});

    runSegemehlAlignment(std::move(args), alignmentsFastqOutPath, "Could not align reads");
}

void SegemehlAligner::alignPairedReads(const fs::path &queryForwardFastqInPath,
                                       const fs::path &queryReverseFastqInPath,
                                       const fs::path &alignmentsFastqOutPath) const {
    const size_t forwardThreads = threadsAdaptedToEntries(queryForwardFastqInPath);
    const size_t reverseThreads = threadsAdaptedToEntries(queryReverseFastqInPath);

    if (forwardThreads == 0 || reverseThreads == 0) {
        writeEmptyAlignments(alignmentsFastqOutPath, forwardThreads == 0 ? queryForwardFastqInPath
                                                                         : queryReverseFastqInPath);
        return;
    }

    const size_t threads = std::min(forwardThreads, reverseThreads);

    auto args = getGeneralAlignmentArgs(threads);

    args.insert(args.end(),
                {"-q", queryForwardFastqInPath.string(), "-p", queryReverseFastqInPath.string(),
                 "-o", alignmentsFastqOutPath.string(), "-H",
                 std::to_string(static_cast<int>(!parameters.multimapAlignments))});

    runSegemehlAlignment(std::move(args), alignmentsFastqOutPath, "Could not align reads");
}

auto SegemehlAligner::convertToCStrings(std::vector<std::string> &args) -> std::vector<char *> {
    std::vector<char *> c_args(args.size() + 1);
    std::ranges::transform(args, c_args.begin(), [](std::string &arg) { return arg.data(); });

    c_args.back() = nullptr;  // argv must be null terminated

    return c_args;
}

}  // namespace pipelines::align
