#include "SingleEndPreprocessor.hpp"

// standard
#include <params/basic.h>

#include <cstddef>
#include <filesystem>
#include <functional>
#include <future>
#include <ranges>
#include <string>
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
#include "PreprocessFilter.hpp"
#include "PreprocessSample.hpp"
#include "RecordTrimmer.hpp"
#include "Utility.hpp"
#include "CheckedFastqReader.hpp"
#include "utility/ConcurrentInput.hpp"

namespace pipelines::preprocess {

using namespace dataTypes;

void SingleEndPreprocessor::process(const PreprocessSampleSingle& sample) const {
    Logger::log("Processing sample (single-end): ", sample.input.sampleName);

    ChunkResult totalResult = parameters.deduplicate ? processWithDeduplication(sample)
                                                     : processWithoutDeduplication(sample);

    Logger::log("Poly-G trimming: ", totalResult.polyGChangedRecords, " records changed, ",
                totalResult.polyGRemovedBases, " bases removed");

    Logger::log("Merging temporary files");

    const auto tmpFilePaths = helper::getFilePathsInDir(sample.output.tmpFastqDir);
    helper::mergeFastqFiles(tmpFilePaths, sample.output.outputFastqPath);
    helper::deleteDir(sample.output.tmpFastqDir);

    Logger::log("Finished processing sample: ", totalResult.getPassedRecords(), " passed records, ",
                totalResult.getFailedRecords(), " failed records");
}

auto SingleEndPreprocessor::processWithDeduplication(const PreprocessSampleSingle& sample) const -> ChunkResult {
    auto retained = Deduplicator::deduplicate(DeduplicationBySequenceSingleConfig{sample.input.inputFastqPath});
    return processInput(sample, &retained.validRecordOrdinals);
}

auto SingleEndPreprocessor::processWithoutDeduplication(const PreprocessSampleSingle& sample) const -> ChunkResult {
    return processInput(sample, nullptr);
}

auto SingleEndPreprocessor::processInput(const PreprocessSampleSingle& sample,
                                       const std::set<size_t>* retained) const -> ChunkResult {
    CheckedFastqReader reader{sample.input.inputFastqPath};
    size_t ordinal = 0;
    utility::ConcurrentInput<FastqRecord> input{parameters.chunkSize, [&]() -> std::optional<FastqRecord> {
        while (auto record = reader.next()) {
            const size_t current = ordinal++;
            if (!retained || retained->contains(current)) return record;
        }
        return std::nullopt;
    }};
    return utility::consumeConcurrently(input, parameters.threadCount - 1,
        [&] { return processChunk(input, sample.output.tmpFastqDir); });
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
            const auto before = record.sequence().size();
            RecordTrimmer::trim3PolyG(record, parameters.minPolyGCount);
            const auto removed = before - record.sequence().size();
            result.polyGChangedRecords += removed != 0;
            result.polyGRemovedBases += removed;
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
    polyGChangedRecords += other.polyGChangedRecords;
    polyGRemovedBases += other.polyGRemovedBases;
    passedRecords += other.passedRecords;
    failedRecords += other.failedRecords;
}

}  // namespace pipelines::preprocess
