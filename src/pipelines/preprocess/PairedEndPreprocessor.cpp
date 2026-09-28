#include "PairedEndPreprocessor.hpp"

// Standard
#include <params/basic.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <future>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>

// Internal
#include "DeduplicationConfig.hpp"
#include "Deduplicator.hpp"
#include "FastqRecord.hpp"
#include "Logger.hpp"
#include "PairedRecordMerger.hpp"
#include "PreprocessFilter.hpp"
#include "PreprocessSample.hpp"
#include "RecordTrimmer.hpp"
#include "Utility.hpp"
#include "CheckedFastqReader.hpp"
#include "utility/ConcurrentInput.hpp"

namespace pipelines::preprocess {

void PairedEndPreprocessor::ChunkResult::operator+=(
    const PairedEndPreprocessor::ChunkResult& other) {
    polyGChangedRecords += other.polyGChangedRecords;
    polyGRemovedBases += other.polyGRemovedBases;
    mergedRecords += other.getMergedRecords();
    singleFwdRecords += other.getSingleFwdRecords();
    singleRevRecords += other.getSingleRevRecords();
    pairedRecords += other.getPairedRecords();
    failedMergedRecords += other.getFailedMergedRecords();
    failedForwardRecords += other.getFailedForwardRecords();
    failedReverseRecords += other.getFailedReverseRecords();
}

void PairedEndPreprocessor::process(const PreprocessSamplePaired& sample) const {
    Logger::log("Processing sample (paired-end): ", sample.input.sampleName);

    ChunkResult totalResult = parameters.deduplicate ? processWithDeduplication(sample)
                                                     : processWithoutDeduplication(sample);

    Logger::log("Poly-G trimming: ", totalResult.polyGChangedRecords, " records changed, ",
                totalResult.polyGRemovedBases, " bases removed");

    Logger::log("Merging temporary files");

    const auto tmpMergedFilePaths = helper::getFilePathsInDir(sample.output.tmpMergedFastqDir);
    helper::mergeFastqFiles(tmpMergedFilePaths, sample.output.outputMergedFastqPath);
    helper::deleteDir(sample.output.tmpMergedFastqDir);

    const auto tmpSingletonForwardFilePaths =
        helper::getFilePathsInDir(sample.output.tmpSingletonForwardFastqDir);
    helper::mergeFastqFiles(tmpSingletonForwardFilePaths,
                            sample.output.outputSingletonForwardFastqPath);
    helper::deleteDir(sample.output.tmpSingletonForwardFastqDir);

    const auto tmpSingletonReverseFilePaths =
        helper::getFilePathsInDir(sample.output.tmpSingletonReverseFastqDir);
    helper::mergeFastqFiles(tmpSingletonReverseFilePaths,
                            sample.output.outputSingletonReverseFastqPath);
    helper::deleteDir(sample.output.tmpSingletonReverseFastqDir);

    auto tmpPairedForwardFilePaths =
        helper::getFilePathsInDir(sample.output.tmpPairedForwardFastqDir);
    std::ranges::sort(tmpPairedForwardFilePaths);
    helper::mergeFastqFiles(tmpPairedForwardFilePaths, sample.output.outputPairedForwardFastqPath);
    helper::deleteDir(sample.output.tmpPairedForwardFastqDir);

    auto tmpPairedReverseFilePaths =
        helper::getFilePathsInDir(sample.output.tmpPairedReverseFastqDir);
    std::ranges::sort(tmpPairedReverseFilePaths);
    helper::mergeFastqFiles(tmpPairedReverseFilePaths, sample.output.outputPairedReverseFastqPath);
    helper::deleteDir(sample.output.tmpPairedReverseFastqDir);

    Logger::log("Finished processing sample: ", totalResult.getMergedRecords(), " merged, ",
                totalResult.getPairedRecords(), " non-merged paired records, ",
                totalResult.getSingleFwdRecords(), " single forward records, ",
                totalResult.getSingleRevRecords(), " single reverse records, ",
                totalResult.getFailedMergedRecords(), " failed merged records, ",
                totalResult.getFailedForwardRecords(), " failed forward records, ",
                totalResult.getFailedReverseRecords(), " failed reverse records");
}

auto PairedEndPreprocessor::processWithDeduplication(const PreprocessSamplePaired& sample) const -> ChunkResult {
    auto retained = Deduplicator::deduplicate(DeduplicationBySequencePairedConfig{sample.input.inputForwardFastqPath, sample.input.inputReverseFastqPath});
    return processInput(sample, &retained.validRecordOrdinals);
}

auto PairedEndPreprocessor::processWithoutDeduplication(const PreprocessSamplePaired& sample) const -> ChunkResult {
    return processInput(sample, nullptr);
}

auto PairedEndPreprocessor::processInput(const PreprocessSamplePaired& sample,
                                       const std::set<size_t>* retained) const -> ChunkResult {
    CheckedFastqPairReader reader{sample.input.inputForwardFastqPath, sample.input.inputReverseFastqPath};
    size_t ordinal = 0;
    utility::ConcurrentInput<PairedFastqRecords> input{parameters.chunkSize / 2, [&]() -> std::optional<PairedFastqRecords> {
        while (auto record = reader.next()) {
            const size_t current = ordinal++;
            if (!retained || retained->contains(current)) return record;
        }
        return std::nullopt;
    }};
    return utility::consumeConcurrently(input, parameters.threadCount - 1,
        [&] { return processChunk(input, sample.output); });
}

void PairedEndPreprocessor::trimWindowedQuality(
    PairedFastqRecords& records, const RecordTrimmer::TrimWindowedConfig& config) const {
    if (parameters.windowTrimmingSize > 0) {
        RecordTrimmer::trimWindowedQuality(records.first, config);
        RecordTrimmer::trimWindowedQuality(records.second, config);
    }
};

void PairedEndPreprocessor::trimAdapters(PairedFastqRecords& records) const {
    for (auto const& adapter : adapters5fwd) {
        RecordTrimmer::trimAdapter(adapter, records.first, parameters.minOverlapTrimming);
    }

    for (auto const& adapter : adapters3fwd) {
        RecordTrimmer::trimAdapter(adapter, records.first, parameters.minOverlapTrimming);
    }

    for (auto const& adapter : adapters5rev) {
        RecordTrimmer::trimAdapter(adapter, records.second, parameters.minOverlapTrimming);
    }

    for (auto const& adapter : adapters3rev) {
        RecordTrimmer::trimAdapter(adapter, records.second, parameters.minOverlapTrimming);
    }
};

template <typename T>
[[nodiscard]] auto PairedEndPreprocessor::processChunk(
    T& recordIterator, const PrepocessSampleOutputPaired& tmpOutDir) const -> ChunkResult {
    ChunkResult result;
    const std::string uuid = helper::getUUID();

    const fs::path tmpMergedFastqOutPath = tmpOutDir.tmpMergedFastqDir / (uuid + ".fastq.gz");
    seqan3::sequence_file_output mergedOut{tmpMergedFastqOutPath};

    const fs::path tmpSingletonFwdFastqOutPath =
        tmpOutDir.tmpSingletonForwardFastqDir / (uuid + ".fastq.gz");
    seqan3::sequence_file_output snglFwdOut{tmpSingletonFwdFastqOutPath};

    const fs::path tmpSingletonRevFastqOutPath =
        tmpOutDir.tmpSingletonReverseFastqDir / (uuid + ".fastq.gz");
    seqan3::sequence_file_output snglRevOut{tmpSingletonRevFastqOutPath};

    const fs::path tmpPairedFwdFastqOutPath =
        tmpOutDir.tmpPairedForwardFastqDir / (uuid + ".fastq.gz");
    seqan3::sequence_file_output pairedFwdOut{tmpPairedFwdFastqOutPath};

    const fs::path tmpPairedRevFastqOutPath =
        tmpOutDir.tmpPairedReverseFastqDir / (uuid + ".fastq.gz");
    seqan3::sequence_file_output pairedRevOut{tmpPairedRevFastqOutPath};

    for (auto&& [record1, record2] : recordIterator) {
        PairedFastqRecords records{std::make_pair(std::move(record1), std::move(record2))};

        if (parameters.trimPolyG) {
            for (auto* record : {&records.first, &records.second}) {
                const auto before = record->sequence().size();
                RecordTrimmer::trim3PolyG(*record, parameters.minPolyGCount);
                const auto removed = before - record->sequence().size();
                result.polyGChangedRecords += removed != 0;
                result.polyGRemovedBases += removed;
            }
        }

        trimWindowedQuality(records, {.windowTrimmingSize = parameters.windowTrimmingSize,
                                      .minMeanWindowPhred = parameters.minMeanWindowQuality});

        trimAdapters(records);

        const PreprocessFilter::Criteria filterCriteria{
            .minLengthThreshold = parameters.minLengthThreshold,
            .minQualityThreshold = parameters.minQualityThreshold};
        const bool filtFwd = PreprocessFilter::passes(filterCriteria, records.first);
        const bool filtRev = PreprocessFilter::passes(filterCriteria, records.second);

        if (filtFwd && filtRev) {
            auto mergedRecord = PairedRecordMerger::mergeRecordPair(
                records, parameters.minOverlapMerging, parameters.maxMissMatchFractionMerging);

            if (!mergedRecord.has_value()) {
                pairedFwdOut.push_back(records.first);
                pairedRevOut.push_back(records.second);

                result.incrementPairedRecords();
                continue;
            }

            if (PreprocessFilter::passes(filterCriteria, mergedRecord.value())) {
                mergedOut.push_back(mergedRecord.value());

                result.incrementMergedRecords();
                continue;
            }

            result.incrementFailedMergedRecords();
            continue;
        }

        if (filtFwd) {
            snglFwdOut.push_back(records.first);

            result.incrementSingleFwdRecords();
        } else {
            result.incrementFailedForwardRecords();
        }

        if (filtRev) {
            snglRevOut.push_back(records.second);

            result.incrementSingleRevRecords();
        } else {
            result.incrementFailedReverseRecords();
        }
    }

    return result;
}

}  // namespace pipelines::preprocess
