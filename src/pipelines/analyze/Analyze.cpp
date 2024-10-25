#include "Analyze.hpp"

// Standard
#include <sys/stat.h>

#include <future>
#include <string>
#include <vector>

// Internal
#include "Constants.hpp"
#include "FeatureWriter.hpp"
#include "InteractionClusterGenerator.hpp"
#include "InteractionsWriter.hpp"
#include "SplitRecordsParser.hpp"
#include "StatisticEvaluator.hpp"

namespace pipelines::analyze {

// Analyze
void Analyze::process(const AnalyzeData &data) {
    Logger::log(LogLevel::INFO, constants::pipelines::PROCESSING_TREATMENT_MESSAGE);

    std::vector<std::future<void>> futures;

    futures.reserve(data.treatmentSamples.size());
    for (const auto &sample : data.treatmentSamples) {
        futures.push_back(std::async(std::launch::async, &Analyze::processSample, this, sample));
    }

    if (!data.controlSamples.has_value()) {
        Logger::log(LogLevel::INFO, "No control samples provided");
        return;
    }

    for (const auto &sample : data.controlSamples.value()) {
        futures.push_back(std::async(std::launch::async, &Analyze::processSample, this, sample));
    }

    for (auto &future : futures) {
        future.get();
    }
}

void Analyze::processSample(AnalyzeSample sample) {
    Logger::log(LogLevel::INFO, "Processing sample: ", sample.input.sampleName);

    std::vector<InteractionCluster> clusters =
        SplitRecordsParser::parse(sample.input.splitAlignmentsPath);

    seqan3::sam_file_input splitsIn{sample.input.splitAlignmentsPath, sam_field_ids{}};

    auto &header = splitsIn.header();
    const std::deque<std::string> &referenceIDs = header.ref_ids();

    InteractionClusterGenerator clusterGenerator{sample.input.sampleName,
                                                 featureAnnotator,
                                                 referenceIDs,
                                                 parameters.featureOrientation,
                                                 parameters.minimumClusterReadCount,
                                                 parameters.clusterDistanceThreshold};

    auto mergingResult = clusterGenerator.mergeClusters(std::move(clusters));

    assignNonAnnotatedContiguousToSupplementaryFeatures(
        sample.input.unassignedContiguousAlignmentsPath,
        mergingResult.supplementaryFeatureAnnotator, mergingResult.featureCounts);

    assignAnnotatedContiguousFragmentCountsToTranscripts(
        sample.input.contiguousAlignmentsTranscriptCountsPath, mergingResult.featureCounts);

    const size_t totalFragmentCount =
        parseSampleFragmentCount(sample.input.sampleFragmentCountsPath);
    const auto evaluatedClusters =
        StatisticEvaluator::evaluate(mergingResult.annotatedClusters, mergingResult.featureCounts,
                                     totalFragmentCount, parameters.padjThreshold);

    writeTranscriptCounts(mergingResult.featureCounts,
                          sample.output.interactionsTranscriptCountsPath);

    annotation::FeatureWriter::write(
        mergingResult.supplementaryFeatureAnnotator.getFeatureTreeMap(),
        sample.output.supplementaryFeaturesPath, annotation::FileType::GFF);

    const InteractionsWriter::OutputPaths outputPaths{
        .interactionsOutputPath = sample.output.interactionsPath,
        .interactionsBEDOutputPath = sample.output.interactionsBEDPath,
        .interactionsBEDArcOutputPath = sample.output.interactionsBEDARCPath};

    InteractionsWriter::writeInteractions(sample.input.sampleName, outputPaths, referenceIDs,
                                          evaluatedClusters);
}

void Analyze::assignAnnotatedContiguousFragmentCountsToTranscripts(
    const fs::path &contiguousTranscriptCountsInPath,
    std::unordered_map<std::string, size_t> &transcriptCounts) {
    std::ifstream transcriptCountsIn(contiguousTranscriptCountsInPath);

    if (!transcriptCountsIn.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", contiguousTranscriptCountsInPath);
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
                transcriptCounts[transcriptID] += std::stoul(token);
            }
            ++column;
        }

        break;
    }
}

void Analyze::assignNonAnnotatedContiguousToSupplementaryFeatures(
    const fs::path &unassignedSingletonsInPath, annotation::FeatureAnnotator &featureAnnotator,
    std::unordered_map<std::string, size_t> &transcriptCounts) {
    seqan3::sam_file_input unassignedSingletonsIn{unassignedSingletonsInPath.string(),
                                                  sam_field_ids{}};

    const auto annotationOrientation = parameters.featureOrientation;

    for (auto &&records : unassignedSingletonsIn) {
        const auto region = dataTypes::GenomicRegion::fromSamRecord(
            records, unassignedSingletonsIn.header().ref_ids());

        if (!region) {
            continue;
        }

        const auto bestFeature =
            featureAnnotator.getBestOverlappingFeature(region.value(), annotationOrientation);

        if (!bestFeature) {
            continue;
        }

        const std::string &transcriptID = bestFeature->id;
        assert(transcriptCounts.contains(transcriptID));
        ++transcriptCounts[transcriptID];
    }
}

auto Analyze::parseSampleFragmentCount(const fs::path &sampleCountsInPath) -> size_t {
    std::ifstream sampleCountsIn(sampleCountsInPath);

    if (!sampleCountsIn.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", sampleCountsInPath.string());
    }

    sampleCountsIn.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // skip header

    size_t totalTranscriptCount = 0;
    std::string line;
    while (std::getline(sampleCountsIn, line)) {
        std::istringstream iss(line);

        size_t column = 0;
        for (std::string token; std::getline(iss, token, '\t');) {
            if (column != 0) {
                totalTranscriptCount += std::stoul(token);
            }
            ++column;
        }

        break;  // Should only be one line
    }

    return totalTranscriptCount;
}

void Analyze::writeTranscriptCounts(const std::unordered_map<std::string, size_t> &featureCounts,
                                    const fs::path &transcriptCountsOutPath) {
    std::ofstream transcriptCountsOut(transcriptCountsOutPath);

    if (!transcriptCountsOut.is_open()) {
        Logger::log(LogLevel::ERROR, "Could not open file: ", transcriptCountsOutPath);
    }

    for (const auto &[transcriptID, count] : featureCounts) {
        transcriptCountsOut << transcriptID << "\t" << count << "\n";
    }
}

}  // namespace pipelines::analyze
