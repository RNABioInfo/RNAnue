#include "Detect.hpp"

#include <cassert>
#include <cstddef>
#include <deque>
#include <filesystem>
#include <functional>
#include <future>
#include <string>
#include <vector>

// Internal
#include "Constants.hpp"
#include "CustomSamTags.hpp"  // NOLINT
#include "DetectSample.hpp"
#include "Logger.hpp"
#include "Utility.hpp"

using namespace dataTypes;

namespace pipelines::detect {

void Detect::process(const DetectData& data) {
    Logger::log(LogLevel::INFO, constants::pipelines::PROCESSING_TREATMENT_MESSAGE);

    for (const auto& sample : data.treatmentSamples) {
        processSample(sample);
    }

    if (!data.controlSamples.has_value()) {
        Logger::log(LogLevel::INFO, "No control samples provided");
        return;
    }

    Logger::log(LogLevel::INFO, constants::pipelines::PROCESSING_CONTROL_MESSAGE);

    for (const auto& sample : data.controlSamples.value()) {
        processSample(sample);
    }
}

void Detect::processSample(const DetectSample& sample) const {
    Logger::log(LogLevel::INFO, "Processing sample: ", sample.input.sampleName);

    const fs::path outputTmpDir = sample.output.outputSplitAlignmentsPath.parent_path() / "tmp";

    const ChunkedOutTmpDirs outTmpDirs = prepareTmpOutputDirs(outputTmpDir);

    Logger::log(LogLevel::DEBUG, "Alignments path: ", sample.input.inputAlignmentsPath);

    seqan3::sam_file_input alignmentsIn{sample.input.inputAlignmentsPath, sam_field_ids{}};

    std::vector<size_t> referenceLengths{};
    std::ranges::transform(alignmentsIn.header().ref_id_info, std::back_inserter(referenceLengths),
                           [](auto const& info) { return std::get<0>(info); });
    const std::deque<std::string>& referenceIDs = alignmentsIn.header().ref_ids();

    assert(referenceLengths.size() == referenceIDs.size());

    AsyncGroupBufferType recordInputBuffer =
        alignmentsIn | AsyncSplitRecordGroupBuffer(params.threadCount + 1);

    std::vector<std::future<Result>> processResults;

    for (size_t i = 1; i < params.threadCount; ++i) {
        processResults.emplace_back(std::async(
            std::launch::async, &Detect::processRecordChunk, this, std::ref(outTmpDirs),
            std::ref(recordInputBuffer), std::ref(referenceIDs), std::ref(referenceLengths)));
    }

    Detect::Result mergedResults;

    for (auto& resultFuture : processResults) {
        mergedResults += resultFuture.get();
    }

    Logger::log(LogLevel::INFO, "Processed ", mergedResults.processedRecordsCount, " reads. Found ",
                mergedResults.splitFragmentsCount, " valid split fragments and ",
                mergedResults.singletonFragmentsCount, " singleton fragments. Removed ",
                mergedResults.removedDueToLowMappingQuality,
                " reads due to low mapping quality and ", mergedResults.removedDueToFragmentLength,
                " reads due to short read length.");

    writeTranscriptCountsFile(sample.output.outputContiguousAlignmentsTranscriptCountsPath,
                              mergedResults.transcriptCounts);

    // Merge results from all processing chunks
    mergeOutputFiles(outTmpDirs, sample.output);

    writeReadCountsSummaryFile(mergedResults, sample.input.sampleName,
                               sample.output.outputSharedReadCountsPath);

    fs::remove_all(outputTmpDir);
}

auto Detect::processRecordChunk(const ChunkedOutTmpDirs& outTmpDirs,
                                AsyncGroupBufferType& recordInputBuffer,
                                const std::deque<std::string>& refIDs,
                                const std::vector<size_t>& refLengths) const -> Detect::Result {
    const std::string chunkID = helper::getUUID();

    const fs::path splitsOutPath = outTmpDirs.outputTmpSplitsDir / (chunkID + ".bam");
    const fs::path multiSplitsOutPath = outTmpDirs.outputTmpMultisplitsDir / (chunkID + ".bam");
    const fs::path unassignedContiguousOutPath =
        outTmpDirs.outputTmpUnassignedContiguousDir / (chunkID + ".bam");

    seqan3::sam_file_output splitsOut{splitsOutPath, refIDs, refLengths, sam_field_ids{}};
    seqan3::sam_file_output multiSplitsOut{multiSplitsOutPath, refIDs, refLengths, sam_field_ids{}};
    seqan3::sam_file_output unassignedContiguousOut{unassignedContiguousOutPath, refIDs, refLengths,
                                                    sam_field_ids{}};

    size_t recordsCount = 0;
    size_t totalSplitFragmentsCount = 0;
    size_t totalSingletonFragmentsCount = 0;

    size_t removedDueToLowMapQuality = 0;
    size_t removedDueToReadLength = 0;

    TranscriptCounts singletonTranscriptCounts;

    auto assignSingletonTranscriptCount = [&](SamRecord& record) {
        if (!record.reference_id() || !record.reference_position()) [[unlikely]] {
            return;
        }

        const auto region = GenomicRegion::fromSamRecord(record, refIDs);

        if (!region) [[unlikely]] {
            return;
        }

        const auto bestFeature =
            featureAnnotator.getBestOverlappingFeature(region.value(), params.featureOrientation);

        if (bestFeature) {
            singletonTranscriptCounts[bestFeature->id]++;
        } else {
            unassignedContiguousOut.push_back(record);
        }
    };

    for (std::vector<SamRecord>& recordGroup : recordInputBuffer) {
        recordsCount += recordGroup.size();

        std::unordered_set<size_t> invalidRecordHitGroups;
        std::vector<SamRecord> validRecordsGroup;

        auto isInvalid = [&params = params](SamRecord& record) {
            return static_cast<bool>(record.flag() & seqan3::sam_flag::unmapped) ||
                   record.mapping_quality() < params.minimumMapQuality ||
                   record.sequence().size() < params.minimumFragmentLength;
        };

        for (SamRecord& record : recordGroup) {
            if (isInvalid(record)) {
                invalidRecordHitGroups.insert(record.tags().get<"HI"_tag>());

                if (record.mapping_quality() < params.minimumMapQuality) {
                    removedDueToLowMapQuality++;
                } else {
                    removedDueToReadLength++;
                }

                continue;
            }

            // Read is singleton and does not contain any split information ->
            // Counted for total fragments
            if (!record.tags().contains("XJ"_tag) || record.tags().get<"XJ"_tag>() < 2) {
                if (static_cast<bool>(record.flag() & seqan3::sam_flag::secondary_alignment)) {
                    continue;
                }

                assignSingletonTranscriptCount(record);
                totalSingletonFragmentsCount++;

                continue;
            }

            // Read is part of a split read and valid but a previous record with the
            // same hit group was invalid -> Skip
            if (invalidRecordHitGroups.contains(record.tags().get<"HI"_tag>())) {
                continue;
            }

            validRecordsGroup.push_back(std::move(record));
        }

        std::erase_if(validRecordsGroup, [&](const SamRecord& record) {
            return invalidRecordHitGroups.contains(record.tags().get<"HI"_tag>());
        });

        size_t fragmentsCount =
            processReadRecords(validRecordsGroup, refIDs, splitsOut, multiSplitsOut);
        totalSplitFragmentsCount += fragmentsCount;
    }

    return {.processedRecordsCount = recordsCount,
            .transcriptCounts = singletonTranscriptCounts,
            .splitFragmentsCount = totalSplitFragmentsCount,
            .singletonFragmentsCount = totalSingletonFragmentsCount,
            .removedDueToLowMappingQuality = removedDueToLowMapQuality,
            .removedDueToFragmentLength = removedDueToReadLength};
}

auto Detect::getSplitRecordsEvaluatorParameters(const DetectParameters& params) const
    -> SplitRecordsEvaluationParameters::ParameterVariant {
    if (params.removeSplicingEvents) {
        return SplitRecordsEvaluationParameters::SplicingParameters{
            .baseParameters = {.minComplementarity = params.minimumComplementarity,
                               .minComplementarityFraction = params.minimumSiteLengthRatio,
                               .mfeThreshold = params.maxHybridizationEnergy},
            .orientation = params.featureOrientation,
            .splicingTolerance = params.splicingTolerance,
            .featureAnnotator = &featureAnnotator};
    }

    return SplitRecordsEvaluationParameters::BaseParameters{
        .minComplementarity = params.minimumComplementarity,
        .minComplementarityFraction = params.minimumSiteLengthRatio,
        .mfeThreshold = params.maxHybridizationEnergy};
}

auto Detect::getReferenceIDs(const fs::path& mappingsInPath) -> std::deque<std::string> {
    seqan3::sam_file_input alignmentsIn{mappingsInPath.string(), sam_field_ids{}};

    return alignmentsIn.header().ref_ids();
}

/**
 * Retrieves the final optimal split records and returns the count of fragments
 * for a specific record id. Expects read records to stem exactly from one
 * record id.
 *
 * @param readRecords The vector of SamRecord objects representing the read
 * records.
 * @param splitsOut The output stream for split records.
 * @param multiSplitsOut The output stream for multi split records.
 * @return The count of fragments for a specific record id.
 */
auto Detect::processReadRecords(const std::vector<SamRecord>& readRecords,
                                const std::deque<std::string>& referenceIDs, auto& splitsOut,
                                auto& multiSplitsOut [[maybe_unused]]) const -> size_t {
    if (readRecords.empty()) {
        return 0;
    }

    const auto splitRecords = getSplitRecords(readRecords, referenceIDs);

    if (!splitRecords.has_value()) {
        return 0;
    }

    // TODO Implement multi split writing
    writeSamFile(splitsOut, splitRecords.value().splitRecords);

    return splitRecords.value().splitRecords.size();
}

auto Detect::constructSplitRecords(const SamRecord& readRecord) const
    -> std::optional<SplitRecords> {
    // Number of expected split records for within the record
    const size_t expectedSplitRecords = readRecord.tags().get<"XH"_tag>();

    SplitRecords splitRecords{};
    splitRecords.reserve(expectedSplitRecords);

    std::vector<seqan3::cigar> currentCigar{};
    constexpr size_t initialCigarCapacity = 10;
    splitRecords.reserve(initialCigarCapacity);

    size_t referencePosition = readRecord.reference_position().value_or(0);
    size_t startPosRead{};
    size_t endPosRead{};  // absolute position in read/alignment (e.g., 1 to end)
    size_t startPosSplit{};
    size_t endPosSplit{};  // position in split read (e.g., XX:i to XY:i / 14 to 20)

    bool isValid = true;
    int nextSplitReferenceShift = 0;

    auto const addSplitRecord = [&]() {
        const auto splitSeq =
            readRecord.sequence() | seqan3::views::slice(static_cast<ptrdiff_t>(startPosRead),
                                                         static_cast<ptrdiff_t>(endPosRead));
        const auto splitQual =
            readRecord.base_qualities() | seqan3::views::slice(static_cast<ptrdiff_t>(startPosRead),
                                                               static_cast<ptrdiff_t>(endPosRead));

        if (splitSeq.size() < params.minimumFragmentLength) {
            isValid = false;
            return;
        }

        seqan3::sam_tag_dictionary tags{};
        tags.get<"XX"_tag>() = static_cast<int>(startPosSplit);
        tags.get<"XY"_tag>() = static_cast<int>(endPosSplit);
        tags.get<"XN"_tag>() = static_cast<float>(splitRecords.size());

        splitRecords.emplace_back(readRecord.id(), readRecord.flag(), readRecord.reference_id(),
                                  referencePosition, readRecord.mapping_quality(), currentCigar,
                                  seqan3::dna5_vector(splitSeq.begin(), splitSeq.end()),
                                  std::vector<seqan3::phred42>(splitQual.begin(), splitQual.end()),
                                  std::move(tags));
    };

    auto const addOtherCigar = [&](const auto& cigar) {
        const auto cigarValue = get<0>(cigar);
        endPosRead += cigarValue;
        endPosSplit += cigarValue;
        currentCigar.push_back(cigar);
    };

    auto const addInsertionCigar = [&](const auto& cigar) {
        const auto cigarValue = get<0>(cigar);
        endPosRead += cigarValue;
        endPosSplit += cigarValue;
        nextSplitReferenceShift -= cigarValue;
        currentCigar.push_back(cigar);
    };

    auto const addDeletionCigar = [&](const auto& cigar) {
        currentCigar.push_back(cigar);
        nextSplitReferenceShift += get<0>(cigar);
    };

    auto const addSoftClipCigar = [&](const auto& cigar) {
        const auto cigarValue = get<0>(cigar);
        if (!params.excludeSoftClipping) {
            nextSplitReferenceShift -= cigarValue;
            addOtherCigar(cigar);
            return;
        }

        /* If current cigar is empty, we are at the beginning of the read in case
        of soft clipping at the end of the read it is just ignored */
        if (currentCigar.empty()) {
            nextSplitReferenceShift -= cigarValue;
            startPosRead += cigarValue;
            endPosRead += cigarValue;
            startPosSplit += cigarValue;
            endPosSplit += cigarValue;
        }
    };

    auto const addSkipCigar = [&](const auto& cigar) {
        if (currentCigar.empty()) {
            return;
        }

        addSplitRecord();

        // Set up positions for the next split
        const auto cigarValue = get<0>(cigar);
        assert((cigarValue + endPosRead + nextSplitReferenceShift) >= 0);
        referencePosition += cigarValue + endPosRead + nextSplitReferenceShift;
        startPosSplit = endPosSplit;
        startPosRead = endPosRead;
        nextSplitReferenceShift = 0;
        currentCigar.clear();
    };

    for (const auto& cigar : readRecord.cigar_sequence()) {
        if (cigar == 'M'_cigar_operation || cigar == '='_cigar_operation ||
            cigar == 'X'_cigar_operation) {
            addOtherCigar(cigar);
        } else if (cigar == 'I'_cigar_operation) {
            addInsertionCigar(cigar);
        } else if (cigar == 'D'_cigar_operation) {
            addDeletionCigar(cigar);
        } else if (cigar == 'S'_cigar_operation) {
            addSoftClipCigar(cigar);
        } else if (cigar == 'N'_cigar_operation) {
            addSkipCigar(cigar);
        }
    }

    addSplitRecord();

    if (!isValid) {
        return std::nullopt;
    }

    if (splitRecords.size() != expectedSplitRecords) {
        Logger::log(LogLevel::WARNING, "Expected ", expectedSplitRecords,
                    " split records, but got ", splitRecords.size(),
                    ". Record ID: ", readRecord.id());

        return std::nullopt;
    }

    return splitRecords;
}

/**
 * Constructs split records for a list of records with one or more elements
 * that comprise a splitted read.
 *
 * @param readRecords The vector of read records.
 * @return An optional containing the split records if construction is
 * successful, otherwise std::nullopt.
 */
auto Detect::constructSplitRecords(const std::vector<SamRecord>& readRecords) const
    -> std::optional<SplitRecords> {
    // Number of expected split records for whole read
    const size_t expectedSplitRecords = readRecords.front().tags().get<"XJ"_tag>();

    SplitRecords splitRecords{};
    splitRecords.reserve(expectedSplitRecords);

    for (const auto& record : readRecords) {
        auto splitRecord = constructSplitRecords(record);

        if (!splitRecord.has_value()) {
            return std::nullopt;
        }

        splitRecords.insert(splitRecords.end(), splitRecord->begin(), splitRecord->end());
    }

    if (splitRecords.size() != expectedSplitRecords) {
        Logger::log(LogLevel::WARNING, "Expected ", expectedSplitRecords,
                    " split records, but got ", splitRecords.size(),
                    ". Record ID: ", readRecords.front().id());

        return std::nullopt;
    }

    return splitRecords;
}

/**
 * Retrieves the split records from the given vector of read records.
 *
 * This function groups the read records based on their "HI" tag and constructs
 * split records for each group. It then evaluates the split records and
 * returns the best evaluated split records.
 *
 * @param readRecords The vector of read records.
 * @return An optional containing the best evaluated split records, or an empty
 * optional if no split records were found.
 */
auto Detect::getSplitRecords(const std::vector<SamRecord>& readRecords,
                             const std::deque<std::string>& referenceIDs) const
    -> std::optional<SplitRecordsEvaluator::EvaluatedSplitRecords> {
    std::unordered_map<size_t, std::vector<SamRecord>> recordHitGroups{};
    for (const auto& record : readRecords) {
        recordHitGroups[record.tags().get<"HI"_tag>()].push_back(record);
    }

    std::optional<SplitRecordsEvaluator::EvaluatedSplitRecords> bestSplitRecords{};

    const auto insertBestSplitRecords = [&](SplitRecords& splitRecords) {
        const auto evaluationResult = splitRecordsEvaluator.evaluate(splitRecords, referenceIDs);

        if (std::holds_alternative<SplitRecordsEvaluator::EvaluatedSplitRecords>(
                evaluationResult)) {
            const auto& evaluatedSplitRecords =
                std::get<SplitRecordsEvaluator::EvaluatedSplitRecords>(evaluationResult);

            if (!bestSplitRecords.has_value() || evaluatedSplitRecords > bestSplitRecords.value()) {
                bestSplitRecords.emplace(
                    std::get<SplitRecordsEvaluator::EvaluatedSplitRecords>(evaluationResult));
            }
        } else {
            Logger::log(LogLevel::DEBUG, "Split records failed evaluation. Reason: ",
                        std::get<SplitRecordsEvaluator::FilterReason>(evaluationResult));
        }
    };

    for (const auto& [_, hitGroup] : recordHitGroups) {
        auto splitRecords = constructSplitRecords(hitGroup);

        if (splitRecords.has_value()) {
            insertBestSplitRecords(splitRecords.value());
        }
    }

    return bestSplitRecords;
}

void Detect::writeSamFile(auto& samOut, const std::vector<SamRecord>& splitRecords) {
    for (auto&& record : splitRecords) {
        auto [id, flag, ref_id, ref_offset, mapq, cigar, seq, qual, tags] = record;
        samOut.emplace_back(id, flag, ref_id, ref_offset, mapq, cigar, seq, qual, tags);
    }
}

auto Detect::prepareTmpOutputDirs(const fs::path& tmpOutDir) -> Detect::ChunkedOutTmpDirs {
    const fs::path outputTmpSplitsDir = tmpOutDir / "tmp_splits";
    const fs::path outputTmpMultisplitsDir = tmpOutDir / "tmp_multisplits";
    const fs::path outputTmpUnassignedContiguousRecordsDir =
        tmpOutDir / "tmp_unassigned_contiguous";

    fs::create_directories(outputTmpSplitsDir);
    fs::create_directories(outputTmpMultisplitsDir);
    fs::create_directories(outputTmpUnassignedContiguousRecordsDir);

    return {outputTmpSplitsDir, outputTmpMultisplitsDir, outputTmpUnassignedContiguousRecordsDir};
}

void Detect::writeReadCountsSummaryFile(const Result& results, const std::string& sampleName,
                                        const fs::path& statsFilePath) {
    std::ofstream statsFileStream(statsFilePath);

    if (!statsFileStream.is_open()) {
        throw std::runtime_error("Could not open the stats file.");
    }

    statsFileStream << "sample\tsplits\tsingletons" << std::endl;
    statsFileStream << sampleName << "\t" << results.splitFragmentsCount << "\t"
                    << results.singletonFragmentsCount << "\n";
}

void Detect::mergeOutputFiles(const ChunkedOutTmpDirs& tmpDirs, const DetectOutput& output) {
    std::vector<fs::path> splitsOutFilePaths =
        helper::getValidFilePaths(tmpDirs.outputTmpSplitsDir, {".bam"});
    std::vector<fs::path> multisplitsOutFilePaths =
        helper::getValidFilePaths(tmpDirs.outputTmpMultisplitsDir, {".bam"});
    std::vector<fs::path> unassignedContiguousRecordsOutFilePaths =
        helper::getValidFilePaths(tmpDirs.outputTmpUnassignedContiguousDir, {".bam"});

    helper::mergeSamFiles(splitsOutFilePaths, output.outputSplitAlignmentsPath);
    helper::mergeSamFiles(multisplitsOutFilePaths, output.outputMultisplitAlignmentsPath);
    helper::mergeSamFiles(unassignedContiguousRecordsOutFilePaths,
                          output.outputUnassignedContiguousAlignmentsPath);
}

void Detect::writeTranscriptCountsFile(const fs::path& transcriptCountsFilePath,
                                       const TranscriptCounts& transcriptCounts) {
    std::ofstream transcriptCountsFileStream(transcriptCountsFilePath);

    if (!transcriptCountsFileStream.is_open()) {
        throw std::runtime_error("Could not open the transcript counts file.");
    }

    for (const auto& [transcript, count] : transcriptCounts) {
        transcriptCountsFileStream << transcript << "\t" << count << "\n";
    }
}

}  // namespace pipelines::detect
