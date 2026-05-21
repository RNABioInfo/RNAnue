#include "InteractionsWriter.hpp"

// Standard
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <format>
#include <fstream>
#include <map>
#include <string>
#include <string_view>
#include <vector>

// Internal
#include "EvaluatedInteractionCluster.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "Utility.hpp"

namespace pipelines::analyze {

void InteractionsWriter::writeInteractions(
    const std::string& sampleName, const OutputPaths& outputPaths,
    const std::deque<std::string>& referenceIDs,
    const std::vector<EvaluatedInteractionCluster>& evaluatedClusters) {
    std::ofstream interactionsOutput(outputPaths.interactionsOutputPath);
    if (!interactionsOutput.is_open()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not open file: ",
                                                            outputPaths.interactionsOutputPath);
    }

    std::ofstream interactionsBEDOutput(outputPaths.interactionsBEDOutputPath);
    if (!interactionsBEDOutput.is_open()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not open file: ",
                                                            outputPaths.interactionsBEDOutputPath);
    }

    std::ofstream interactionsBEDArcOutput(outputPaths.interactionsBEDArcOutputPath);
    if (!interactionsBEDArcOutput.is_open()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
            "Could not open file: ", outputPaths.interactionsBEDArcOutputPath);
    }

    std::ofstream interactionReadIDsOutput(outputPaths.interactionReadIDsOutputPath);
    if (!interactionReadIDsOutput.is_open()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
            "Could not open file: ", outputPaths.interactionReadIDsOutputPath);
    }

    writeInteractionsHeader(interactionsOutput);
    writeInteractionsBEDHeader(interactionsBEDOutput, sampleName);
    writeInteractionsBEDArcHeader(interactionsBEDArcOutput, sampleName);
    writeInteractionReadIDsHeader(interactionReadIDsOutput);
    writeInteractionArmCoverageBedGraph(sampleName, evaluatedClusters, referenceIDs,
                                        outputPaths.interactionArmCoverageBedGraphOutputPath);

    size_t intramolecularCount = 0;
    size_t intermolecularCount = 0;

    size_t clusterID = 0;
    for (const auto& cluster : evaluatedClusters) {
        if (cluster.getFirstFeatureID() == cluster.getSecondFeatureID()) {
            ++intramolecularCount;
        } else {
            ++intermolecularCount;
        }

        const std::string clusterIDString = std::to_string(clusterID);

        writeInteraction(cluster, clusterIDString, referenceIDs, interactionsOutput);
        writeInteractionBED(cluster, clusterIDString, referenceIDs, interactionsBEDOutput);
        writeInteractionBEDArc(cluster, clusterIDString, referenceIDs, interactionsBEDArcOutput);
        writeInteractionReadIDs(cluster, clusterIDString, interactionReadIDsOutput);

        ++clusterID;
    }

    Logger::log("After filtering kept ", intramolecularCount + intermolecularCount,
                " split interactions");
    Logger::log("Of which ", intramolecularCount, " are intramolecular interactions");
    Logger::log("Of which ", intermolecularCount, " are intermolecular interactions");
}

void InteractionsWriter::writeInteractionsHeader(std::ofstream& interactionsOut) {
    interactionsOut
        << "cluster_ID\tfst_feat_id\tfst_seg_chr\tfst_seg_strt\tfst_seg_"
           "end\tfst_seg_strd\tsec_feat_id\t"
           "sec_seg_chr\tsec_seg_strt\tsec_seg_end\tsec_seg_strd\tno_splits\t"
           "total_span_bp\teffective_coverage_span_bp\tsupport_per_total_bp\t"
           "support_per_effective_bp\tcoverage_concentration\tcoverage_components\tarm_balance\t"
           "coverage_profile\tmean_inter_crosslinks\tsd_inter_crosslinks\t"
           "gcs\tghs\tp_value\tpadj_value\n";
}

void InteractionsWriter::writeInteractionsBEDHeader(std::ofstream& interactionsBEDOut,
                                                    const std::string& sampleName) {
    interactionsBEDOut << "track name=\"" << sampleName
                       << " RNA-RNA interactions\" description=\"Segments of interacting RNA "
                          "clusters derived from DDD-Experiment\" itemRgb=\"On\"\n";
}

void InteractionsWriter::writeInteractionsBEDArcHeader(std::ofstream& interactionsBEDArcOut,
                                                       const std::string& sampleName) {
    interactionsBEDArcOut
        << "track graphType=arc name=\"" << sampleName
        << " RNA-RNA interactions arcs\" description=\"Segments of interacting RNA "
           "clusters derived from DDD-Experiment\" itemRgb=\"On\"\n";
}

void InteractionsWriter::writeInteractionReadIDsHeader(std::ofstream& interactionReadIDsOut) {
    interactionReadIDsOut << "cluster_ID\tread_IDs\n";
}

void InteractionsWriter::writeInteraction(const EvaluatedInteractionCluster& cluster,
                                          const std::string& clusterID,
                                          const std::deque<std::string>& referenceIDs,
                                          std::ofstream& interactionOut) {
    interactionOut << "cluster" << clusterID << "\t";

    interactionOut << cluster.getFirstFeatureID() << "\t";
    interactionOut << getReferenceID(cluster.getFirstSegment().getReferenceIDIndex(), referenceIDs)
                   << "\t";
    interactionOut << cluster.getFirstSegment().getStart() << "\t";
    interactionOut << cluster.getFirstSegment().getEnd() << "\t";
    interactionOut << static_cast<char>(cluster.getFirstSegment().getStrand()) << "\t";

    interactionOut << cluster.getSecondFeatureID() << "\t";
    interactionOut << getReferenceID(cluster.getSecondSegment().getReferenceIDIndex(), referenceIDs)
                   << "\t";
    interactionOut << cluster.getSecondSegment().getStart() << "\t";
    interactionOut << cluster.getSecondSegment().getEnd() << "\t";
    interactionOut << static_cast<char>(cluster.getSecondSegment().getStrand()) << "\t";

    const CoverageShapeMetrics coverageMetrics = cluster.coverageShapeMetrics();

    interactionOut << std::format("{:.2f}", cluster.getTranscriptContribution()) << "\t";
    interactionOut << coverageMetrics.totalSpanBp << "\t";
    interactionOut << std::format("{:.2f}", coverageMetrics.effectiveCoverageSpanBp) << "\t";
    interactionOut << std::format("{:.6f}", coverageMetrics.supportPerTotalBp) << "\t";
    interactionOut << std::format("{:.6f}", coverageMetrics.supportPerEffectiveBp) << "\t";
    interactionOut << std::format("{:.4f}", coverageMetrics.coverageConcentration) << "\t";
    interactionOut << coverageMetrics.coverageComponents << "\t";
    interactionOut << std::format("{:.4f}", coverageMetrics.armBalance) << "\t";
    interactionOut << coverageMetrics.coverageProfile << "\t";
    interactionOut << std::format("{:.2f}", cluster.meanCrosslinkingSiteCount()) << "\t";
    interactionOut << std::format("{:.2f}", cluster.standardDeviationCrosslinkingSiteCount())
                   << "\t";
    interactionOut << std::format("{:.2f}", cluster.complementarityStatistics()) << "\t";
    interactionOut << std::format("{:.2f}", cluster.hybridizationEnergyStatistics()) << "\t";
    interactionOut << std::format("{:.4f}", cluster.getPValue()) << "\t";
    interactionOut << std::format("{:.4f}", cluster.getPadj());
    interactionOut << "\n";
}

void InteractionsWriter::writeInteractionBED(const EvaluatedInteractionCluster& cluster,
                                             const std::string& clusterID,
                                             const std::deque<std::string>& referenceIDs,
                                             std::ofstream& bedOut) {
    const std::string color = helper::generateRandomHexColor();

    bedOut << getReferenceID(cluster.getFirstSegment().getReferenceIDIndex(), referenceIDs) << "\t";
    bedOut << cluster.getFirstSegment().getStart() << "\t";
    bedOut << cluster.getFirstSegment().getEnd() << "\t";
    bedOut << "cluster" << clusterID << "_segment1"
           << "\t";
    bedOut << "0\t";
    bedOut << static_cast<char>(cluster.getFirstSegment().getStrand()) << "\t";
    bedOut << cluster.getFirstSegment().getStart() << "\t";
    bedOut << cluster.getFirstSegment().getEnd() << "\t";
    bedOut << color << "\n";

    bedOut << getReferenceID(cluster.getSecondSegment().getReferenceIDIndex(), referenceIDs)
           << "\t";
    bedOut << cluster.getSecondSegment().getStart() << "\t";
    bedOut << cluster.getSecondSegment().getEnd() << "\t";
    bedOut << "cluster" << clusterID << "_segment2"
           << "\t";
    bedOut << "0\t";
    bedOut << static_cast<char>(cluster.getSecondSegment().getStrand()) << "\t";
    bedOut << cluster.getSecondSegment().getStart() << "\t";
    bedOut << cluster.getSecondSegment().getEnd() << "\t";
    bedOut << color << "\n";
}

void InteractionsWriter::writeInteractionBEDArc(const EvaluatedInteractionCluster& cluster,
                                                const std::string& clusterID,
                                                const std::deque<std::string>& referenceIDs,
                                                std::ofstream& bedArcOut) {
    bedArcOut << getReferenceID(cluster.getFirstSegment().getReferenceIDIndex(), referenceIDs)
              << "\t";
    bedArcOut << cluster.getFirstSegment().getStart() << "\t";
    bedArcOut << cluster.getSecondSegment().getStart() << "\t";
    bedArcOut << "cluster" << clusterID << "\n";
}

void InteractionsWriter::writeInteractionArmCoverageBedGraph(
    const std::string& sampleName,
    const std::vector<EvaluatedInteractionCluster>& evaluatedClusters,
    const std::deque<std::string>& referenceIDs,
    const fs::path& interactionArmCoverageBedGraphOutputPath) {
    struct CoverageEvent {
        int32_t position{};
        double delta{};
    };

    std::map<int32_t, std::vector<CoverageEvent>> coverageEventsByReference;
    auto addCoverageRuns = [&coverageEventsByReference](int32_t referenceIDIndex,
                                                        const std::vector<CoverageRun>& runs) {
        auto& events = coverageEventsByReference[referenceIDIndex];
        events.reserve(events.size() + runs.size() * 2);

        for (const auto& run : runs) {
            if (run.start >= run.end || run.coverage <= 0.0) {
                continue;
            }

            events.push_back({.position = run.start, .delta = run.coverage});
            events.push_back({.position = run.end, .delta = -run.coverage});
        }
    };

    for (const auto& cluster : evaluatedClusters) {
        addCoverageRuns(cluster.getFirstSegment().getReferenceIDIndex(),
                        cluster.getFirstArmCoverageRuns());
        addCoverageRuns(cluster.getSecondSegment().getReferenceIDIndex(),
                        cluster.getSecondArmCoverageRuns());
    }

    std::ofstream bedGraphOut(interactionArmCoverageBedGraphOutputPath);
    if (!bedGraphOut.is_open()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
            "Could not open file: ", interactionArmCoverageBedGraphOutputPath);
    }

    bedGraphOut << "track type=bedGraph name=\"" << sampleName
                << " weighted interaction arm coverage\" description=\"Weighted per-base "
                   "coverage across retained RNAnue interaction arms\"\n";

    constexpr double coverageEpsilon = 1e-12;
    for (const auto& [referenceIDIndex, events] : coverageEventsByReference) {
        std::vector<CoverageEvent> sortedEvents = events;
        std::ranges::sort(sortedEvents, [](const CoverageEvent& lhs, const CoverageEvent& rhs) {
            return lhs.position < rhs.position;
        });

        double currentCoverage = 0.0;
        int32_t previousPosition = sortedEvents.empty() ? 0 : sortedEvents.front().position;
        size_t index = 0;

        while (index < sortedEvents.size()) {
            const int32_t position = sortedEvents[index].position;

            if (position > previousPosition && currentCoverage > coverageEpsilon) {
                bedGraphOut << getReferenceID(referenceIDIndex, referenceIDs) << "\t"
                            << previousPosition << "\t" << position << "\t"
                            << std::format("{:.6f}", currentCoverage) << "\n";
            }

            double delta = 0.0;
            while (index < sortedEvents.size() && sortedEvents[index].position == position) {
                delta += sortedEvents[index].delta;
                ++index;
            }

            currentCoverage += delta;
            if (std::abs(currentCoverage) <= coverageEpsilon) {
                currentCoverage = 0.0;
            }
            previousPosition = position;
        }
    }
}

void InteractionsWriter::writeInteractionReadIDs(const EvaluatedInteractionCluster& cluster,
                                                 const std::string& clusterID,
                                                 std::ofstream& interactionReadIDsOut) {
    interactionReadIDsOut << "cluster" << clusterID << "\t";
    constexpr std::string_view seperator = ",";

    std::string_view curSeperator{};
    for (const auto& readID : cluster.getRecordIDs()) {
        interactionReadIDsOut << std::format("{}{}", curSeperator, readID);
        curSeperator = seperator;
    }

    interactionReadIDsOut << "\n";
}
}  // namespace pipelines::analyze
