#include "PairedEndPreprocessor.hpp"

// Standard
#include <params/basic.h>

#include <future>
#include <utility>

// Internal
#include "DeduplicationConfig.hpp"
#include "Deduplicator.hpp"
#include "FastqRecord.hpp"
#include "PairedRecordMerger.hpp"
#include "PreprocessFilter.hpp"
#include "RecordTrimmer.hpp"
#include "Utility.hpp"

namespace pipelines::preprocess {

void PairedEndPreprocessor::ChunkResult::operator+=(
    const PairedEndPreprocessor::ChunkResult& other) {
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

[[nodiscard]] auto PairedEndPreprocessor::processWithDeduplication(
    const PreprocessSamplePaired& sample) const -> ChunkResult {
    std::vector<std::future<ChunkResult>> processResults;
    processResults.reserve(parameters.threadCount - 1);

    const auto deduplicationConfig =
        DeduplicationBySequencePairedConfig{.recordsPathFwd = sample.input.inputForwardFastqPath,
                                            .recordsPathRev = sample.input.inputReverseFastqPath};
    auto deduplicatedResults = Deduplicator::deduplicate(deduplicationConfig);

    const size_t totalRecords = deduplicatedResults.recordPairs.size();
    const size_t numThreads = parameters.threadCount;
    const size_t chunkSize = (totalRecords + numThreads - 1) / numThreads;  // Ceiling division

    std::vector<std::vector<std::pair<FastqRecord, FastqRecord>>> chunks;
    chunks.reserve(numThreads);

    for (size_t i = 0; i < numThreads; ++i) {
        size_t start = i * chunkSize;
        size_t end = std::min(start + chunkSize, totalRecords);
        if (start >= end) {
            break;
        }

        std::vector<std::pair<FastqRecord, FastqRecord>> chunk{
            std::make_move_iterator(deduplicatedResults.recordPairs.begin() +
                                    static_cast<long>(start)),
            std::make_move_iterator(deduplicatedResults.recordPairs.begin() +
                                    static_cast<long>(end))};

        chunks.push_back(std::move(chunk));
    }

    const size_t numChunks = chunks.size();

    for (size_t i = 1; i < numChunks; ++i) {
        processResults.emplace_back(std::async(
            std::launch::async,
            &PairedEndPreprocessor::processChunk<std::vector<std::pair<FastqRecord, FastqRecord>>>,
            this, std::ref(chunks[i]), std::cref(sample.output)));
    }

    ChunkResult totalResult = processChunk(chunks[0], sample.output);

    for (auto& resultFuture : processResults) {
        totalResult += resultFuture.get();
    }

    return totalResult;
};

[[nodiscard]] auto PairedEndPreprocessor::processWithoutDeduplication(
    const PreprocessSamplePaired& sample) const -> ChunkResult {
    seqan3::sequence_file_input recForwardIn{sample.input.inputForwardFastqPath};
    seqan3::sequence_file_input recReverseIn{sample.input.inputReverseFastqPath};

    PairedEndAsyncInputBuffer pairedInputBuffer =
        seqan3::views::zip(recForwardIn, recReverseIn) |
        seqan3::views::async_input_buffer(parameters.chunkSize / 2);

    std::vector<std::future<ChunkResult>> processResults;
    processResults.reserve(parameters.threadCount - 1);

    for (size_t i = 1; i < parameters.threadCount; ++i) {
        processResults.emplace_back(std::async(
            std::launch::async, &PairedEndPreprocessor::processChunk<PairedEndAsyncInputBuffer>,
            this, std::ref(pairedInputBuffer), std::cref(sample.output)));
    }

    ChunkResult totalResult = {};

    for (auto& result : processResults) {
        totalResult += result.get();
    }

    return totalResult;
};

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
            RecordTrimmer::trim3PolyG(record1);
            RecordTrimmer::trim3PolyG(record2);
        }

        trimWindowedQuality(records, {.windowTrimmingSize = parameters.windowTrimmingSize,
                                      .minMeanWindowPhred = parameters.minMeanWindowQuality});

        trimAdapters(records);

        const PreprocessFilter::Criteria filterCriteria{
            .minLengthThreshold = parameters.minLengthThreshold,
            .minQualityThreshold = parameters.minQualityThreshold};
        const bool filtFwd = PreprocessFilter::passes(filterCriteria, record1);
        const bool filtRev = PreprocessFilter::passes(filterCriteria, record2);

        if (filtFwd && filtRev) {
            auto mergedRecord = PairedRecordMerger::mergeRecordPair(
                records, parameters.minOverlapMerging, parameters.maxMissMatchFractionMerging);

            if (!mergedRecord.has_value()) {
                pairedFwdOut.push_back(std::move(records.first));
                pairedRevOut.push_back(std::move(records.second));

                result.incrementPairedRecords();
                continue;
            }

            if (PreprocessFilter::passes(filterCriteria, (mergedRecord.value()))) {
                mergedOut.push_back(std::move(mergedRecord.value()));

                result.incrementMergedRecords();
                continue;
            }

            result.incrementFailedMergedRecords();
            continue;
        }

        if (filtFwd) {
            snglFwdOut.push_back(std::move(records.first));

            result.incrementSingleFwdRecords();
        } else {
            result.incrementFailedForwardRecords();
        }

        if (filtRev) {
            snglRevOut.push_back(std::move(records.second));

            result.incrementSingleRevRecords();
        } else {
            result.incrementFailedReverseRecords();
        }
    }

    return result;
}

}  // namespace pipelines::preprocess
