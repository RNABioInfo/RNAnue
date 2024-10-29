#include "Preprocess.hpp"

#include <future>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include "PairedRecordMerger.hpp"
#include "RecordTrimmer.hpp"
#include "Utility.hpp"

namespace pipelines::preprocess {
Preprocess::Preprocess(PreprocessParameters params) : parameters(std::move(params)) {}

void Preprocess::process(const PreprocessData &data) const {
    Logger::log(LogLevel::INFO, constants::pipelines::PROCESSING_TREATMENT_MESSAGE);

    for (const auto &sample : data.treatmentSamples) {
        processSample(sample);
    }

    if (!data.controlSamples.has_value()) {
        Logger::log(LogLevel::INFO, "No control samples provided");
        return;
    }

    Logger::log(LogLevel::INFO, constants::pipelines::PROCESSING_CONTROL_MESSAGE);

    for (const auto &sample : data.controlSamples.value()) {
        processSample(sample);
    }
}

void Preprocess::processSample(const PreprocessSampleType &sample) const {
    std::visit(
        overloaded{[this](const PreprocessSampleSingle &sample) { processSingleEnd(sample); },
                   [this](const PreprocessSamplePaired &sample) { processPairedEnd(sample); }},
        sample);
}

/**
 * Processes a single-end sample.
 *
 * @param sample The sample to be processed.
 */
void Preprocess::processSingleEnd(const PreprocessSampleSingle &sample) const {
    Logger::log(LogLevel::INFO, "Processing sample (single-end): ", sample.input.sampleName);

    const std::vector<Adapter> adapters5 =
        Adapter::loadAdapters(parameters.adapter5Forward, parameters.maxMissMatchFractionTrimming,
                              TrimConfig::Mode::FIVE_PRIME);
    const std::vector<Adapter> adapters3 =
        Adapter::loadAdapters(parameters.adapter3Forward, parameters.maxMissMatchFractionTrimming,
                              TrimConfig::Mode::THREE_PRIME);

    seqan3::sequence_file_input recIn{sample.input.inputFastqPath};
    SingleEndAsyncInputBuffer asyncInputBuffer =
        recIn | seqan3::views::async_input_buffer(parameters.chunkSize);

    std::vector<std::future<SingleEndResult>> processResults;

    for (size_t i = 1; i < parameters.threadCount; ++i) {
        processResults.emplace_back(
            std::async(std::launch::async, &Preprocess::processSingleEndRecordChunk, this,
                       std::ref(asyncInputBuffer), std::ref(adapters5), std::ref(adapters3),
                       std::ref(sample.output.tmpFastqDir)));
    }

    SingleEndResult totalResult;

    for (auto &result : processResults) {
        const auto res = result.get();
        totalResult += res;
    }

    const auto tmpFilePaths = helper::getFilePathsInDir(sample.output.tmpFastqDir);
    helper::mergeFastqFiles(tmpFilePaths, sample.output.outputFastqPath);
    helper::deleteDir(sample.output.tmpFastqDir);

    Logger::log(LogLevel::INFO, "Finished processing sample: ", sample.input.sampleName, " (",
                totalResult.passedRecords, " passed, ", totalResult.failedRecords, " failed)");
}

/**
 * Process a paired-end sample.
 *
 * @param sample The sample to be processed.
 */
void Preprocess::processPairedEnd(const PreprocessSamplePaired &sample) const {
    Logger::log(LogLevel::INFO, "Processing sample (paired-end): ", sample.input.sampleName);

    std::vector<Adapter> adapters5f =
        Adapter::loadAdapters(parameters.adapter5Forward, parameters.maxMissMatchFractionTrimming,
                              TrimConfig::Mode::FIVE_PRIME);
    std::vector<Adapter> adapters3f =
        Adapter::loadAdapters(parameters.adapter3Forward, parameters.maxMissMatchFractionTrimming,
                              TrimConfig::Mode::THREE_PRIME);
    std::vector<Adapter> adapters5r =
        Adapter::loadAdapters(parameters.adapter5Reverse, parameters.maxMissMatchFractionTrimming,
                              TrimConfig::Mode::FIVE_PRIME);
    std::vector<Adapter> adapters3r =
        Adapter::loadAdapters(parameters.adapter3Reverse, parameters.maxMissMatchFractionTrimming,
                              TrimConfig::Mode::THREE_PRIME);

    seqan3::sequence_file_input recForwardIn{sample.input.inputForwardFastqPath};
    seqan3::sequence_file_input recReverseIn{sample.input.inputReverseFastqPath};

    PairedEndAsyncInputBuffer pairedInputBuffer =
        seqan3::views::zip(recForwardIn, recReverseIn) |
        seqan3::views::async_input_buffer(parameters.chunkSize / 2);

    std::vector<std::future<PairedEndResult>> processResults;

    for (size_t i = 1; i < parameters.threadCount; ++i) {
        processResults.emplace_back(
            std::async(std::launch::async, &Preprocess::processPairedEndRecordChunk, this,
                       std::ref(pairedInputBuffer), std::ref(adapters5f), std::ref(adapters3f),
                       std::ref(adapters5r), std::ref(adapters3r), std::ref(sample.output)));
    }

    PairedEndResult totalResult;

    for (auto &result : processResults) {
        const auto res = result.get();
        totalResult += res;
    }

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

    Logger::log(LogLevel::INFO, "Finished processing sample: ", sample.input.sampleName, " (",
                totalResult.mergedRecords, " merged, ", totalResult.singleFwdRecords,
                " single forward, ", totalResult.singleRevRecords, " single reverse,\n",
                totalResult.failedMergedRecords, " failed merged, ",
                totalResult.failedForwardRecords, " failed forward, ",
                totalResult.failedReverseRecords, " failed reverse)");
}

/**
 * Checks if a record passes the quality and length filters.
 *
 * @param record The record to be checked.
 * @return True if the record passes all the filters, false otherwise.
 */
auto Preprocess::passesFilters(const auto &record) const -> bool {
    // Filter for mean quality
    const auto phredQual =
        record.base_qualities() | std::views::transform([](auto qual) { return qual.to_phred(); });
    const double sum = std::accumulate(phredQual.begin(), phredQual.end(), 0);
    const double meanPhred = sum / std::ranges::size(phredQual);
    const bool passesQual = meanPhred >= double(parameters.minQualityThreshold);

    // Filter for length
    const bool passesLen = record.sequence().size() >= parameters.minLengthThreshold;

    return passesQual && passesLen;
}

auto Preprocess::processSingleEndRecordChunk(SingleEndAsyncInputBuffer &asyncInputBuffer,
                                             const std::vector<Adapter> &adapters5,
                                             const std::vector<Adapter> &adapters3,
                                             const fs::path &tmpOutDir) const
    -> Preprocess::SingleEndResult {
    const std::string uuid = helper::getUUID();
    fs::path tmpFastqOutPath = tmpOutDir / (uuid + ".fastq.gz");

    seqan3::sequence_file_output recOut{tmpFastqOutPath};

    size_t failedRecords = 0;
    size_t passedRecords = 0;

    for (auto &record : asyncInputBuffer) {
        if (parameters.trimPolyG) {
            RecordTrimmer::trim3PolyG(record);
        }

        if (parameters.windowTrimmingSize > 0) {
            RecordTrimmer::trimWindowedQuality(record, parameters.windowTrimmingSize,
                                               parameters.minMeanWindowQuality);
        }

        for (auto const &adapter : adapters5) {
            RecordTrimmer::trimAdapter(adapter, record, parameters.minOverlapTrimming);
        }

        for (auto const &adapter : adapters3) {
            RecordTrimmer::trimAdapter(adapter, record, parameters.minOverlapTrimming);
        }

        if (!passesFilters(record)) {
            ++failedRecords;
            continue;
        }

        recOut.push_back(record);
        ++passedRecords;
    }

    return {.passedRecords = passedRecords, .failedRecords = failedRecords};
}

Preprocess::PairedEndResult Preprocess::processPairedEndRecordChunk(
    Preprocess::PairedEndAsyncInputBuffer &pairedRecordInputBuffer,
    const std::vector<Adapter> &adapters5f, const std::vector<Adapter> &adapters3f,
    const std::vector<Adapter> &adapters5r, const std::vector<Adapter> &adapters3r,
    const PrepocessSampleOutputPaired &sampleOutput) const {
    PairedEndResult result;
    const std::string uuid = helper::getUUID();

    fs::path tmpMergedFastqOutPath = sampleOutput.tmpMergedFastqDir / (uuid + ".fastq.gz");
    seqan3::sequence_file_output mergedOut{tmpMergedFastqOutPath};

    fs::path tmpSingletonFwdFastqOutPath =
        sampleOutput.tmpSingletonForwardFastqDir / (uuid + ".fastq.gz");
    seqan3::sequence_file_output snglFwdOut{tmpSingletonFwdFastqOutPath};

    fs::path tmpSingletonRevFastqOutPath =
        sampleOutput.tmpSingletonReverseFastqDir / (uuid + ".fastq.gz");
    seqan3::sequence_file_output snglRevOut{tmpSingletonRevFastqOutPath};

    for (auto &&[record1, record2] : pairedRecordInputBuffer) {
        if (parameters.trimPolyG) {
            RecordTrimmer::trim3PolyG(record1);
            RecordTrimmer::trim3PolyG(record2);
        }

        if (parameters.windowTrimmingSize > 0) {
            RecordTrimmer::trimWindowedQuality(record1, parameters.windowTrimmingSize,
                                               parameters.minMeanWindowQuality);
            RecordTrimmer::trimWindowedQuality(record2, parameters.windowTrimmingSize,
                                               parameters.minMeanWindowQuality);
        }

        for (auto const &adapter : adapters5f) {
            RecordTrimmer::trimAdapter(adapter, record1, parameters.minOverlapTrimming);
        }

        for (auto const &adapter : adapters3f) {
            RecordTrimmer::trimAdapter(adapter, record1, parameters.minOverlapTrimming);
        }

        for (auto const &adapter : adapters5r) {
            RecordTrimmer::trimAdapter(adapter, record2, parameters.minOverlapTrimming);
        }

        for (auto const &adapter : adapters3r) {
            RecordTrimmer::trimAdapter(adapter, record2, parameters.minOverlapTrimming);
        }

        const bool filtFwd = passesFilters(record1);
        const bool filtRev = passesFilters(record2);

        if (filtFwd && filtRev) {
            auto mergedRecord =
                PairedRecordMerger::mergeRecordPair(record1, record2, parameters.minOverlapMerging,
                                                    parameters.maxMissMatchFractionMerging);

            if (!mergedRecord.has_value()) {
                snglFwdOut.push_back(std::move(record1));
                snglRevOut.push_back(std::move(record2));

                result.singleFwdRecords += 1;
                result.singleRevRecords += 1;

                continue;
            }

            const auto &mergedRecordValue = mergedRecord.value();

            if (passesFilters(mergedRecordValue)) {
                mergedOut.push_back(std::move(mergedRecord.value()));

                result.mergedRecords += 1;
                continue;
            }

            result.failedMergedRecords += 1;

            continue;
        }

        if (filtFwd) {
            snglFwdOut.push_back(std::move(record1));

            result.singleFwdRecords += 1;
        } else {
            result.failedForwardRecords += 1;
        }

        if (filtRev) {
            snglRevOut.push_back(std::move(record2));
        } else {
            result.failedReverseRecords += 1;
        }
    }

    return result;
}

}  // namespace pipelines::preprocess
