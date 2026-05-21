#include "Analyze.hpp"

// Standard
#include <sys/stat.h>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <deque>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <limits>
#include <memory>
#include <optional>
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
#include "AnnotationFilePicker.hpp"
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

    const auto annotationFilePath = utility::AnnotationFilePicker::getFile(parameters);

    const ReferenceIDToIndexMap referenceIDToIndex =
        annotation::loadReferenceIDToIndexMap(data.getInputFilePaths());
    std::shared_ptr<const FeatureAnnotator> featureAnnotator{
        std::make_shared<const FeatureAnnotator>(annotationFilePath, referenceIDToIndex,
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

    TranscriptContributionsByID backgroundContributions;

    parseAnnotatedContiguousFragmentCountsToTranscripts(
        sample.input.contiguousAlignmentsTranscriptCountsPath, backgroundContributions);

    parseNonAnnotatedContiguousToSupplementaryFeatures(
        sample.input.unassignedContiguousAlignmentsPath,
        mergingResult.supplementaryFeatureAnnotator, backgroundContributions);

    auto evaluatedClusters =
        StatisticEvaluator::evaluate(mergingResult.annotatedClusters, backgroundContributions,
                                     parameters.padjThreshold);
    evaluatedClusters = filterByCoverageMetrics(std::move(evaluatedClusters));

    writeTranscriptCounts(backgroundContributions, sample.output.interactionsTranscriptCountsPath);

    annotation::FeatureWriter::write(
        mergingResult.supplementaryFeatureAnnotator.getFeatureTreeMap(), referenceIDs,
        sample.output.supplementaryFeaturesPath, FileType::GFF);

    const InteractionsWriter::OutputPaths outputPaths{
        .interactionsOutputPath = sample.output.interactionsPath,
        .interactionReadIDsOutputPath = sample.output.interactionsReadIDsPath,
        .interactionsBEDOutputPath = sample.output.interactionsBEDPath,
        .interactionsBEDArcOutputPath = sample.output.interactionsBEDARCPath,
        .interactionArmCoverageBedGraphOutputPath =
            sample.output.interactionsArmCoverageBedGraphPath};

    InteractionsWriter::writeInteractions(sample.input.sampleName, outputPaths, referenceIDs,
                                          evaluatedClusters);
}

auto Analyze::filterByCoverageMetrics(std::vector<EvaluatedInteractionCluster> &&clusters) const
    -> std::vector<EvaluatedInteractionCluster> {
    if (!parameters.minimumSupportPerEffectiveBp && !parameters.maximumCoverageComponents &&
        !parameters.minimumArmBalance) {
        return std::move(clusters);
    }

    std::vector<EvaluatedInteractionCluster> filteredClusters;
    filteredClusters.reserve(clusters.size());

    size_t supportDensityFilteredCount = 0;
    size_t coverageComponentsFilteredCount = 0;
    size_t armBalanceFilteredCount = 0;

    for (auto &cluster : clusters) {
        const CoverageShapeMetrics metrics = cluster.coverageShapeMetrics();
        bool keepCluster = true;

        if (parameters.minimumSupportPerEffectiveBp &&
            metrics.supportPerEffectiveBp < *parameters.minimumSupportPerEffectiveBp) {
            keepCluster = false;
            ++supportDensityFilteredCount;
        }

        if (parameters.maximumCoverageComponents &&
            metrics.coverageComponents > *parameters.maximumCoverageComponents) {
            keepCluster = false;
            ++coverageComponentsFilteredCount;
        }

        if (parameters.minimumArmBalance && metrics.armBalance < *parameters.minimumArmBalance) {
            keepCluster = false;
            ++armBalanceFilteredCount;
        }

        if (keepCluster) {
            filteredClusters.emplace_back(std::move(cluster));
        }
    }

    Logger::log("Coverage-shaped support filters kept ", filteredClusters.size(), " of ",
                clusters.size(), " interactions");
    if (parameters.minimumSupportPerEffectiveBp) {
        Logger::log("Filtered ", supportDensityFilteredCount,
                    " interactions below minimum support_per_effective_bp: ",
                    *parameters.minimumSupportPerEffectiveBp);
    }
    if (parameters.maximumCoverageComponents) {
        Logger::log("Filtered ", coverageComponentsFilteredCount,
                    " interactions above maximum coverage_components: ",
                    *parameters.maximumCoverageComponents);
    }
    if (parameters.minimumArmBalance) {
        Logger::log("Filtered ", armBalanceFilteredCount,
                    " interactions below minimum arm_balance: ", *parameters.minimumArmBalance);
    }

    return filteredClusters;
}

void Analyze::parseAnnotatedContiguousFragmentCountsToTranscripts(
    const fs::path &contiguousTranscriptCountsInPath,
    TranscriptContributionsByID &transcriptCounts) {
    std::ifstream transcriptCountsIn(contiguousTranscriptCountsInPath);

    if (!transcriptCountsIn.is_open()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not open file: ",
                                                            contiguousTranscriptCountsInPath);
    }

    parseTranscriptContributionStream(transcriptCountsIn, transcriptCounts);
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
        const float currentContribution = transcriptCounts[transcriptID];

        if (std::isnan(currentContribution) || std::isnan(contributionScore)) {
            Logger::log<LogLevel::WARNING>(std::format(
                "Invalid contribution – check unassigned contiguous record in detect step. Record "
                "ID: {}, "
                "Interaction ID: {}, Current "
                "contribution: {}, Record contribution: {}",
                record.id(), transcriptID, currentContribution, contributionScore));

            continue;
        }

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
