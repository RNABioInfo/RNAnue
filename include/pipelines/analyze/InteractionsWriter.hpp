#pragma once

// Standard
#include <deque>
#include <filesystem>
#include <string>

// Internal
#include "EvaluatedInteractionCluster.hpp"

namespace pipelines::analyze {

namespace fs = std::filesystem;

class InteractionsWriter {
   public:
    InteractionsWriter() = delete;

    struct OutputPaths;

    static void writeInteractions(
        const std::string& sampleName, const OutputPaths& outputPaths,
        const std::deque<std::string>& referenceIDs,
        const std::vector<EvaluatedInteractionCluster>& evaluatedClusters);

   private:
    [[nodiscard]] static auto getReferenceID(const int32_t referenceIDIndex,
                                             const std::deque<std::string>& referenceIDs)
        -> std::string {
        return referenceIDs[referenceIDIndex];
    }

    static void writeInteractionsHeader(std::ofstream& interactionsOutputFile);
    static void writeInteractionsBEDHeader(std::ofstream& interactionsBEDOut,
                                           const std::string& sampleName);
    static void writeInteractionsBEDArcHeader(std::ofstream& interactionsBEDArcOut,
                                              const std::string& sampleName);

    static void writeInteraction(const EvaluatedInteractionCluster& cluster,
                                 const std::string& clusterID,
                                 const std::deque<std::string>& referenceIDs,
                                 std::ofstream& interactionOut);
    static void writeInteractionBED(const EvaluatedInteractionCluster& cluster,
                                    const std::string& clusterID,
                                    const std::deque<std::string>& referenceIDs,
                                    std::ofstream& bedOut);
    static void writeInteractionBEDArc(const EvaluatedInteractionCluster& cluster,
                                       const std::string& clusterID,
                                       const std::deque<std::string>& referenceIDs,
                                       std::ofstream& bedArcOut);
};

struct InteractionsWriter::OutputPaths {
    fs::path interactionsOutputPath;
    fs::path interactionsBEDOutputPath;
    fs::path interactionsBEDArcOutputPath;
};

}  // namespace pipelines::analyze
