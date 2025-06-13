#include "Analyze.hpp"

// Standard
#include <sys/stat.h>

#include <cassert>
#include <cstddef>
#include <deque>
#include <filesystem>
#include <fstream>
#include <ios>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

// Internal
#include "AnalyzeData.hpp"
#include "AnalyzeParameters.hpp"
#include "AnalyzeSample.hpp"
#include "Constants.hpp"
#include "CustomSamTags.hpp"  // IWYU pragma: keep
#include "FeatureAnnotator.hpp"
#include "FeatureWriter.hpp"
#include "FileType.hpp"
#include "InteractionCluster.hpp"
#include "InteractionsWriter.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "ParallelInteractionClusterGenerator.hpp"
#include "SamRecord.hpp"
#include "SplitRecordsParser.hpp"
#include "StatisticEvaluator.hpp"
#include "TranscriptContributionsByID.hpp"

namespace pipelines::analyze {

using namespace annotation;

// Analyze
void Analyze::process(const AnalyzeData &data) {
    Logger::log(constants::pipelines::PROCESSING_TREATMENT_MESSAGE);

    const ReferenceIDToIndexMap referenceIDToIndex =
        annotation::loadReferenceIDToIndexMap(data.getInputFilePaths());
    std::shared_ptr<const FeatureAnnotator> featureAnnotator{
        std::make_shared<const FeatureAnnotator>(parameters.featuresInPath, referenceIDToIndex,
                                                 parameters.featureTypes)};

    for (const auto &sample : data.treatmentSamples) {
        processSample(sample, featureAnnotator);
    }

    if (!data.controlSamples.has_value()) {
        Logger::log("No control samples provided");
        return;
    }

    for (const auto &sample : data.controlSamples.value()) {
        processSample(sample, featureAnnotator);
    }
}

void Analyze::processSample(const AnalyzeSample &sample,
                            std::shared_ptr<const FeatureAnnotator> featureAnnotator) {
    Logger::log("Processing sample: ", sample.input.sampleName);

    std::vector<InteractionCluster> clusters =
        SplitRecordsParser::parse(sample.input.splitAlignmentsPath);

    seqan3::sam_file_input splitsIn{sample.input.splitAlignmentsPath, SamFieldIDs{}};

    auto &header = splitsIn.header();
    const std::deque<std::string> &referenceIDs = header.ref_ids();

    ParallelInteractionClusterGenerator clusterGenerator{std::move(featureAnnotator),
                                                         parameters.getClusteringParameters()};

    auto mergingResult = clusterGenerator.mergeClusters(std::move(clusters), parameters.threadCount,
                                                        parameters.chunkSize);

    parseNonAnnotatedContiguousToSupplementaryFeatures(
        sample.input.unassignedContiguousAlignmentsPath,
        mergingResult.supplementaryFeatureAnnotator, mergingResult.featureCounts);

    parseAnnotatedContiguousFragmentCountsToTranscripts(
        sample.input.contiguousAlignmentsTranscriptCountsPath, mergingResult.featureCounts);

    const float totalFragmentCount =
        parseSampleFragmentCount(sample.input.sampleFragmentCountsPath);
    const auto evaluatedClusters =
        StatisticEvaluator::evaluate(mergingResult.annotatedClusters, mergingResult.featureCounts,
                                     totalFragmentCount, parameters.padjThreshold);

    writeTranscriptCounts(mergingResult.featureCounts,
                          sample.output.interactionsTranscriptCountsPath);

    annotation::FeatureWriter::write(
        mergingResult.supplementaryFeatureAnnotator.getFeatureTreeMap(), referenceIDs,
        sample.output.supplementaryFeaturesPath, FileType::GFF);

    const InteractionsWriter::OutputPaths outputPaths{
        .interactionsOutputPath = sample.output.interactionsPath,
        .interactionReadIDsOutputPath = sample.output.interactionsReadIDsPath,
        .interactionsBEDOutputPath = sample.output.interactionsBEDPath,
        .interactionsBEDArcOutputPath = sample.output.interactionsBEDARCPath};

    InteractionsWriter::writeInteractions(sample.input.sampleName, outputPaths, referenceIDs,
                                          evaluatedClusters);
}

void Analyze::parseAnnotatedContiguousFragmentCountsToTranscripts(
    const fs::path &contiguousTranscriptCountsInPath,
    TranscriptContributionsByID &transcriptCounts) {
    std::ifstream transcriptCountsIn(contiguousTranscriptCountsInPath);

    if (!transcriptCountsIn.is_open()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not open file: ",
                                                            contiguousTranscriptCountsInPath);
    }

    std::string line;
    while (std::getline(transcriptCountsIn, line)) {
        std::istringstream iss(line);

        size_t column = 0;
        std::string transcriptID;
        for (std::string token; std::getline(iss, token, '\t');) {
            if (column == 0) {
                transcriptID = token;
            } else if (column == 1) {
                transcriptCounts[transcriptID] += std::stof(token);
            }
            ++column;
        }

        break;
    }
}

void Analyze::parseNonAnnotatedContiguousToSupplementaryFeatures(
    const fs::path &unassignedSingletonsInPath,
    const annotation::FeatureAnnotator &featureAnnotator,
    TranscriptContributionsByID &transcriptCounts) {
    seqan3::sam_file_input unassignedSingletonsIn{unassignedSingletonsInPath.string(),
                                                  SamFieldIDs{}};

    const auto annotationOrientation = parameters.featureOrientation;

    for (const auto &record : unassignedSingletonsIn) {
        const auto bestFeature =
            featureAnnotator.getBestOverlappingFeature(record, annotationOrientation);

        if (!bestFeature) {
            continue;
        }

        const std::string &transcriptID = bestFeature->getAnnotationID();
        const float contributionScore = record.tags().get<"XB"_tag>();

        assert(transcriptCounts.contains(transcriptID));
        transcriptCounts[transcriptID] += contributionScore;
    }
}

auto Analyze::parseSampleFragmentCount(const fs::path &sampleCountsInPath) -> float {
    std::ifstream sampleCountsIn(sampleCountsInPath);

    if (!sampleCountsIn.is_open()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not open file: ",
                                                            sampleCountsInPath.string());
    }

    sampleCountsIn.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // skip header

    float totalTranscriptCount = 0;
    std::string line;
    while (std::getline(sampleCountsIn, line)) {
        std::istringstream iss(line);

        size_t column = 0;
        for (std::string token; std::getline(iss, token, '\t');) {
            if (column != 0) {  // Skip name column
                totalTranscriptCount += std::stof(token);
            }
            ++column;
        }

        break;  // Should only be one line
    }

    return totalTranscriptCount;
}

void Analyze::writeTranscriptCounts(const TranscriptContributionsByID &featureCounts,
                                    const fs::path &transcriptCountsOutPath) {
    std::ofstream transcriptCountsOut(transcriptCountsOutPath);

    if (!transcriptCountsOut.is_open()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not open file: ",
                                                            transcriptCountsOutPath);
    }

    for (const auto &[transcriptID, count] : featureCounts) {
        transcriptCountsOut << transcriptID << "\t" << count << "\n";
    }
}

}  // namespace pipelines::analyze
