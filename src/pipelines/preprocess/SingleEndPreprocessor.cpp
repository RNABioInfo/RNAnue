#include "SingleEndPreprocessor.hpp"

// standard
#include <params/basic.h>

#include <future>
#include <string>

// Internal
#include "DeduplicationConfig.hpp"
#include "Deduplicator.hpp"
#include "FastqRecord.hpp"
#include "PreprocessFilter.hpp"
#include "PreprocessSample.hpp"
#include "RecordTrimmer.hpp"

namespace pipelines::preprocess {

using namespace dataTypes;

void SingleEndPreprocessor::process(const PreprocessSampleSingle& sample) const {
    Logger::log("Processing sample (single-end): ", sample.input.sampleName);

    ChunkResult totalResult = parameters.deduplicate ? processWithDeduplication(sample)
                                                     : processWithoutDeduplication(sample);

    Logger::log("Merging temporary files");

    const auto tmpFilePaths = helper::getFilePathsInDir(sample.output.tmpFastqDir);
    helper::mergeFastqFiles(tmpFilePaths, sample.output.outputFastqPath);
    helper::deleteDir(sample.output.tmpFastqDir);

    Logger::log("Finished processing sample: ", totalResult.getPassedRecords(), " passed records, ",
                totalResult.getFailedRecords(), " failed records");
}

auto SingleEndPreprocessor::processWithDeduplication(const PreprocessSampleSingle& sample) const
    -> ChunkResult {
    const auto deduplicationConfig =
        DeduplicationBySequenceSingleConfig{sample.input.inputFastqPath};
    auto deduplicatedResults = Deduplicator::deduplicate(deduplicationConfig);

    auto isValidRecord = [&deduplicatedResults](const FastqRecord& record) {
        return deduplicatedResults.validRecordIDs.contains(record.id());
    };

    std::vector<std::future<ChunkResult>> processResults;
    processResults.reserve(parameters.threadCount - 1);

    seqan3::sequence_file_input recIn{sample.input.inputFastqPath};
    auto asyncInputBuffer = recIn | seqan3::views::async_input_buffer(parameters.chunkSize) |
                            std::ranges::views::filter(isValidRecord);

    for (size_t i = 1; i < parameters.threadCount; ++i) {
        processResults.emplace_back(std::async(
            std::launch::async, &SingleEndPreprocessor::processChunk<decltype(asyncInputBuffer)>,
            this, std::ref(asyncInputBuffer), std::cref(sample.output.tmpFastqDir)));
    }

    // Process data in main thread
    ChunkResult totalResult = processChunk(asyncInputBuffer, sample.output.tmpFastqDir);

    // Collect results from other threads
    for (auto& resultFuture : processResults) {
        totalResult += resultFuture.get();
    }

    return totalResult;
}

auto SingleEndPreprocessor::processWithoutDeduplication(const PreprocessSampleSingle& sample) const
    -> ChunkResult {
    std::vector<std::future<ChunkResult>> processResults;
    processResults.reserve(parameters.threadCount - 1);

    seqan3::sequence_file_input recIn{sample.input.inputFastqPath};
    SingleEndAsyncInputBuffer asyncInputBuffer =
        recIn | seqan3::views::async_input_buffer(parameters.chunkSize);

    for (size_t i = 1; i < parameters.threadCount; ++i) {
        processResults.emplace_back(std::async(
            std::launch::async, &SingleEndPreprocessor::processChunk<SingleEndAsyncInputBuffer>,
            this, std::ref(asyncInputBuffer), std::cref(sample.output.tmpFastqDir)));
    }

    // Process data in main thread
    ChunkResult totalResult = processChunk(asyncInputBuffer, sample.output.tmpFastqDir);

    // Collect results from other threads
    for (auto& resultFuture : processResults) {
        totalResult += resultFuture.get();
    }

    return totalResult;
}

template <typename T>
auto SingleEndPreprocessor::processChunk(T& recordIterator, const fs::path& tmpOutDir) const
    -> ChunkResult {
    ChunkResult result{};

    const std::string uuid = helper::getUUID();
    fs::path tmpFastqOutPath = tmpOutDir / (uuid + ".fastq.gz");

    seqan3::sequence_file_output recOut{tmpFastqOutPath};

    for (auto& record : recordIterator) {
        if (parameters.trimPolyG) {
            RecordTrimmer::trim3PolyG(record);
        }

        if (parameters.windowTrimmingSize > 0) {
            RecordTrimmer::trimWindowedQuality(
                record, {.windowTrimmingSize = parameters.windowTrimmingSize,
                         .minMeanWindowPhred = parameters.minMeanWindowQuality});
        }

        for (auto const& adapter : adapters5) {
            RecordTrimmer::trimAdapter(adapter, record, parameters.minOverlapTrimming);
        }

        for (auto const& adapter : adapters3) {
            RecordTrimmer::trimAdapter(adapter, record, parameters.minOverlapTrimming);
        }

        if (!PreprocessFilter::passes({.minLengthThreshold = parameters.minLengthThreshold,
                                       .minQualityThreshold = parameters.minQualityThreshold},
                                      record)) {
            result.incrementFailedRecords();
            continue;
        }

        recOut.push_back(record);
        result.incrementPassedRecords();
    }

    return result;
}

void SingleEndPreprocessor::ChunkResult::operator+=(const ChunkResult& other) {
    passedRecords += other.passedRecords;
    failedRecords += other.failedRecords;
}

}  // namespace pipelines::preprocess
