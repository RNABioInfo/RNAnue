#include "InteractionsWriter.hpp"

// Standard
#include <cstddef>
#include <deque>
#include <fstream>
#include <string>
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

    writeInteractionsHeader(interactionsOutput);
    writeInteractionsBEDHeader(interactionsBEDOutput, sampleName);
    writeInteractionsBEDArcHeader(interactionsBEDArcOutput, sampleName);

    size_t intramolecularCount = 0;
    size_t intermolecularCount = 0;

    size_t clusterID = 0;
    for (const auto& cluster : evaluatedClusters) {
        if (cluster.getFirstFeatureID() == cluster.getSecondFeatureID()) {
            ++intramolecularCount;
        } else {
            ++intermolecularCount;
        }

        writeInteraction(cluster, std::to_string(clusterID), referenceIDs, interactionsOutput);
        writeInteractionBED(cluster, std::to_string(clusterID), referenceIDs,
                            interactionsBEDOutput);
        writeInteractionBEDArc(cluster, std::to_string(clusterID), referenceIDs,
                               interactionsBEDArcOutput);

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
           "sec_seg_chr\tsec_seg_strt\tsec_seg_end\tsec_seg_strd\tno_splits\tmean_inter_"
           "crosslinks\tsd_inter_crosslinks\t"
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

    interactionOut << cluster.fragmentCount() << "\t";
    interactionOut << cluster.meanCrosslinkingSiteCount() << "\t";
    interactionOut << cluster.standardDeviationCrosslinkingSiteCount() << "\t";
    interactionOut << cluster.complementarityStatistics() << "\t";
    interactionOut << cluster.hybridizationEnergyStatistics() << "\t";
    interactionOut << cluster.getPValue() << "\t";
    interactionOut << cluster.getPadj();
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

}  // namespace pipelines::analyze
