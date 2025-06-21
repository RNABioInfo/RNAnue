#include "Detect.hpp"

#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <future>
#include <ios>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

// seqan3
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

// Internal
#include "AnnotationStep.hpp"
#include "AsyncSplitReadGroupBuffer.hpp"
#include "ComplementarityEvaluationStep.hpp"
#include "Constants.hpp"
#include "DetectData.hpp"
#include "DetectParameters.hpp"
#include "DetectSample.hpp"
#include "EvaluationContextResultHandler.hpp"
#include "FeatureAnnotator.hpp"
#include "HitGroupEvaluator.hpp"
#include "HybridizationEvaluationStep.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "MetricsTrackingStep.hpp"
#include "ReadGroupEvaluationParameters.hpp"
#include "ReadGroupEvaluator.hpp"
#include "ReadGroupPostScoringStep.hpp"
#include "SamRecord.hpp"
#include "SamReference.hpp"
#include "SegemehlReadGroupPreprocessor.hpp"
#include "SplicingEvaluationStep.hpp"
#include "SplitRecords.hpp"
#include "TempOutputDirs.hpp"
#include "Utility.hpp"
#include "VariantOverload.hpp"
#include "seqan3/utility/tuple/pod_tuple.hpp"

using namespace dataTypes;

namespace pipelines::detect {

void Detect::process(const DetectData& data) {
    Logger::log(constants::pipelines::PROCESSING_TREATMENT_MESSAGE);

    const ReferenceIDToIndexMap referenceIDToIndex =
        annotation::loadReferenceIDToIndexMap(data.getInputFilePaths());
    std::shared_ptr<const FeatureAnnotator> featureAnnotator{
        std::make_shared<const FeatureAnnotator>(params.featuresInPath, referenceIDToIndex,
                                                 params.featureTypes)};

    const auto evaluationParameters =
        ReadGroupEvaluationParameters::makeParams(params, featureAnnotator);

    auto processSamples = [&](const auto& evaluationParams) {
        for (const auto& sample : data.treatmentSamples) {
            processSample(sample, evaluationParams);
        }

        if (!data.controlSamples.has_value()) {
            Logger::log("No control samples provided");
            return;
        }

        Logger::log(constants::pipelines::PROCESSING_CONTROL_MESSAGE);

        for (const auto& sample : data.controlSamples.value()) {
            processSample(sample, evaluationParams);
        }
    };

    std::visit(overloaded{[&](const ReadGroupEvaluationParameters::Base& baseParams) {
                              processSamples(baseParams);
                          },
                          [&](const ReadGroupEvaluationParameters::Splicing& splicingParams) {
                              processSamples(splicingParams);
                          }},
               evaluationParameters);
}

template <ReadGroupEvaluationParameters::Type ParamT>
void Detect::processSample(const DetectSample& sample, const ParamT& evaluationParams) const {
    Logger::log("Processing sample: ", sample.input.sampleName);

    // Prepare temporary dirs for chunked processing
    const fs::path outputParentDir = sample.output.outputSplitAlignmentsPath.parent_path();
    const fs::path outputTmpDir = outputParentDir / "tmp";
    const TempOutputDirs outTmpDirs = prepareTmpOutputDirs(outputTmpDir);

    // Get sam file input infos and buffer
    seqan3::sam_file_input alignmentsIn{sample.input.inputAlignmentsPath, SamFieldIDs{}};
    SamReference reference{alignmentsIn.header()};
    AsyncGroupBufferType recordInputBuffer =
        alignmentsIn | AsyncSplitReadGroupBuffer(params.threadCount + 1);

    // Process and store results in a chunked process
    std::vector<std::future<Result>> results;

    for (size_t i = 1; i < params.threadCount; ++i) {
        results.emplace_back(std::async(std::launch::async,
                                        &Detect::template processRecordChunk<ParamT>, this,
                                        std::cref(outTmpDirs), std::ref(recordInputBuffer),
                                        std::ref(reference), std::cref(evaluationParams)));
    }

    Detect::Result mergedResults =
        processRecordChunk(outTmpDirs, recordInputBuffer, reference, evaluationParams);

    for (auto& resultFuture : results) {
        mergedResults += resultFuture.get();
    }

    mergedResults.createPlots(outputParentDir);

    Logger::log(mergedResults.toString());

    writeTranscriptCountsFile(sample.output.outputContiguousAlignmentsTranscriptCountsPath,
                              mergedResults.singletonTranscriptCounts);

    // Merge results from all processing chunks
    mergeTmpFiles(outTmpDirs, sample.output, reference);

    writeReadCountsSummaryFile(mergedResults, sample.input.sampleName,
                               sample.output.outputSharedReadCountsPath);

    fs::remove_all(outputTmpDir);
}

template <ReadGroupEvaluationParameters::Type ParamT>
auto Detect::processRecordChunk(const TempOutputDirs& outTmpDirs,
                                AsyncGroupBufferType& recordInputBuffer, SamReference& reference,
                                const ParamT& evaluationParams) const -> Detect::Result {
    EvaluationContextResultHandler resultHandler{outTmpDirs, reference, params.chunkSize};

    SegemehlReadGroupPreprocessorConfig collationConfig{
        .constructionParams = {.minimumMapQuality = params.minimumMapQuality,
                               .minimumFragmentLength = params.minimumFragmentLength,
                               .excludeSoftClipping = params.excludeSoftClipping},
        .maxPrimaryAlignmentCount = params.maxPrimaryAlignmentCount};
    SegemehlReadGroupPreprocessor preprocessor{collationConfig};

    ReadGroupPostScoringStep postprocessor{{static_cast<float>(params.minHitGroupContribution)}};

    auto complementarityStep = ComplementarityEvaluationStep{
        ComplementarityEvaluationStepConfig::makeConfig(evaluationParams)};

    auto hybridizationStep = HybridizationEvaluationStep{
        HybridizationEvaluationStepConfig::makeConfig(evaluationParams)};

    auto annotationStep = AnnotationStep{AnnotationStepConfig::makeConfig(evaluationParams)};

    auto metricTrackingStep = MetricsTrackingStep{};

    if constexpr (std::same_as<ParamT, ReadGroupEvaluationParameters::Base>) {
        auto evaluator = ReadGroupEvaluator(
            std::ref(preprocessor), std::ref(postprocessor), std::ref(complementarityStep),
            std::ref(hybridizationStep), std::ref(annotationStep), std::ref(metricTrackingStep));

        for (ReadGroup readGroup : recordInputBuffer) {
            auto res = evaluator.evaluate(std::move(readGroup));

            for (const auto& context : res.successes) {
                std::visit(resultHandler, context);
            }
        }

        Logger::log("Finished batch with metrics: ", metricTrackingStep);
    } else {
        auto splicingStep =
            SplicingEvaluationStep{SplicingEvaluationStepConfig::makeConfig(evaluationParams)};

        auto evaluator = ReadGroupEvaluator(std::ref(preprocessor), std::ref(postprocessor),
                                            std::ref(splicingStep), std::ref(complementarityStep),
                                            std::ref(hybridizationStep), std::ref(annotationStep),
                                            std::ref(metricTrackingStep));

        for (ReadGroup readGroup : recordInputBuffer) {
            auto res = evaluator.evaluate(std::move(readGroup));

            for (const auto& context : res.successes) {
                std::visit(resultHandler, context);
            }
        }
    }

    resultHandler.save();

    return {.preprocessMetrics = preprocessor.getMetrics(),
            .complementarityMetrics = complementarityStep.getMetrics(),
            .hybridizationMetrics = hybridizationStep.getMetrics(),
            .postprocessMetrics = postprocessor.getMetrics(),
            .singletonTranscriptCounts = resultHandler.getSingletonTranscriptCounts()};
}

auto Detect::getReferenceIDs(const fs::path& mappingsInPath) -> std::deque<std::string> {
    seqan3::sam_file_input alignmentsIn{mappingsInPath.string(), SamFieldIDs{}};

    return alignmentsIn.header().ref_ids();
}

auto Detect::prepareTmpOutputDirs(const fs::path& tmpOutDir) -> TempOutputDirs {
    const fs::path outputTmpSplitsDir = tmpOutDir / "tmp_splits";
    const fs::path outputTmpMultisplitsDir = tmpOutDir / "tmp_multisplits";
    const fs::path outputTmpUnassignedContiguousRecordsDir =
        tmpOutDir / "tmp_unassigned_contiguous";

    fs::create_directories(outputTmpSplitsDir);
    fs::create_directories(outputTmpMultisplitsDir);
    fs::create_directories(outputTmpUnassignedContiguousRecordsDir);

    return {.outputTmpSplitsDir = outputTmpSplitsDir,
            .outputTmpMultisplitsDir = outputTmpMultisplitsDir,
            .outputTmpUnassignedSingletonDir = outputTmpUnassignedContiguousRecordsDir};
}

void Detect::writeReadCountsSummaryFile(const Result& results, const std::string& sampleName,
                                        const fs::path& statsFilePath) {
    std::ofstream statsFileStream(statsFilePath);

    if (!statsFileStream.is_open()) {
        throw std::runtime_error("Could not open the stats file.");
    }

    const auto& contributionScoreByRecordType =
        results.postprocessMetrics.getContributionScoreByRecordType();

    if (!contributionScoreByRecordType.contains(SplitRecordType::SINGLETON) ||
        !contributionScoreByRecordType.contains(SplitRecordType::CHIMERIC)) {
        Logger::log<SourceLocation{}, LogLevel::ERROR>(
            "Could not retrieve contribution scores for hit group types.");
    }

    statsFileStream << "sample\tsplits\tsingletons\n";
    statsFileStream << sampleName << "\t" << std::fixed
                    << contributionScoreByRecordType.at(SplitRecordType::CHIMERIC) << "\t"
                    << contributionScoreByRecordType.at(SplitRecordType::SINGLETON) << "\n";
}

void Detect::mergeTmpFiles(const TempOutputDirs& tmpDirs, const DetectOutput& output,
                           const SamReference& reference) {
    std::vector<fs::path> splitsOutFilePaths =
        helper::getValidFilePaths(tmpDirs.outputTmpSplitsDir, {".bam"});
    std::vector<fs::path> multisplitsOutFilePaths =
        helper::getValidFilePaths(tmpDirs.outputTmpMultisplitsDir, {".bam"});
    std::vector<fs::path> unassignedContiguousRecordsOutFilePaths =
        helper::getValidFilePaths(tmpDirs.outputTmpUnassignedSingletonDir, {".bam"});

    helper::mergeSamFiles(splitsOutFilePaths, output.outputSplitAlignmentsPath, reference);
    helper::mergeSamFiles(multisplitsOutFilePaths, output.outputMultisplitAlignmentsPath,
                          reference);
    helper::mergeSamFiles(unassignedContiguousRecordsOutFilePaths,
                          output.outputUnassignedContiguousAlignmentsPath, reference);
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
