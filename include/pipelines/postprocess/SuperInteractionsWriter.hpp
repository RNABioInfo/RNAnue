#pragma once

// Standard
#include <algorithm>
#include <cstddef>
#include <format>
#include <fstream>
#include <iterator>
#include <ostream>
#include <ranges>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Internal
#include "Interaction.hpp"
#include "InteractionGenerator.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "PostprocessData.hpp"
#include "PostprocessSample.hpp"
#include "Utility.hpp"
#include "csv.hpp"

namespace pipelines::postprocess::SuperInteractionsWriter {

namespace {
[[nodiscard]] auto getReferenceID(const int referenceIDIndex,
                                  const std::unordered_map<int, std::string>& referenceIDs)
    -> std::string {
    return referenceIDs.at(referenceIDIndex);
}

[[nodiscard]] auto getSampleIDs(
    const std::unordered_map<std::string, std::vector<PostprocessSample>>& samplesMap)
    -> std::vector<std::string> {
    std::vector<std::string> sampleIDs;

    for (const auto& [key, sampleGroup] : samplesMap) {
        for (const auto& sample : sampleGroup) {
            sampleIDs.emplace_back(sample.input.sampleName);
        }
    }

    return sampleIDs;
}

void writeSuperInteractionsGCTHeader(const std::vector<std::string>& sampleIDs,
                                     size_t superInteractionCount, std::ostream& out,
                                     csv::TSVWriter<std::ofstream>& writer) {
    out << std::format("#1.3\n{}\t{}\n", superInteractionCount, sampleIDs.size());
    std::vector<std::string> headerTokens = {"Name", "Description"};
    std::ranges::copy(sampleIDs, std::back_inserter(headerTokens));

    writer << headerTokens;
}

void writeSuperInteractionsBEDPEHeader(std::ofstream& bedOut,
                                       const std::vector<std::string>& sampleNames) {
    bedOut << "#track name=\"Super RNA-RNA interactions of samples: "
           << helper::toString(sampleNames)
           << "\" description=\"Segments of interacting RNA "
              "clusters derived from DDD-Experiment\" itemRgb=\"On\"\n";
    bedOut << "#columns color=12\n";
}

void writeSuperInteractionGCT(csv::TSVWriter<std::ofstream>& writer,
                              const std::string& interactionID,
                              const std::vector<std::string>& sampleIDs,
                              const Interaction& superInteraction) {
    std::vector<std::string> lineTokens{interactionID, interactionID};
    lineTokens.reserve(sampleIDs.size() + 1);

    for (const auto& sampleName : sampleIDs) {
        lineTokens.emplace_back(
            std::format("{:.2f}", superInteraction.getContributionScore(sampleName)));
    }

    writer << lineTokens;
}

void writeSuperInteractionsBEDPE(csv::TSVWriter<std::ofstream>& writer,
                                 const std::unordered_map<int, std::string>& referenceIndexToIDMap,
                                 const std::string& interactionID,
                                 const std::vector<std::string>& sampleIDs,
                                 const Interaction& superInteraction) {
    std::vector<std::string> lineTokens;
    static constexpr int bedpeFieldCount = 12;
    lineTokens.reserve(bedpeFieldCount);

    lineTokens.emplace_back(getReferenceID(superInteraction.getFirstSegment().getReferenceIDIndex(),
                                           referenceIndexToIDMap));
    lineTokens.emplace_back(std::format("{}", superInteraction.getFirstSegment().getStart()));
    lineTokens.emplace_back(std::format("{}", superInteraction.getFirstSegment().getEnd()));
    lineTokens.emplace_back(getReferenceID(
        superInteraction.getSecondSegment().getReferenceIDIndex(), referenceIndexToIDMap));
    lineTokens.emplace_back(std::format("{}", superInteraction.getSecondSegment().getStart()));
    lineTokens.emplace_back(std::format("{}", superInteraction.getSecondSegment().getEnd()));
    lineTokens.emplace_back(interactionID);
    lineTokens.emplace_back(".");
    lineTokens.emplace_back(
        std::string{static_cast<char>(superInteraction.getFirstSegment().getStrand())});
    lineTokens.emplace_back(
        std::string{static_cast<char>(superInteraction.getSecondSegment().getStrand())});

    std::string parentClusters;
    for (const auto& sampleID : sampleIDs) {
        auto clusterIdsView =
            superInteraction.getInteractionIDs() |
            std::views::filter(
                [&sampleID](const auto& interaction) { return interaction.sampleID == sampleID; }) |
            std::views::transform(
                [](const auto& interaction) -> std::string_view { return interaction.clusterID; });

        std::string parentIDs;
        bool first = true;
        for (const auto& clusterID : clusterIdsView) {
            if (!first) {
                parentIDs.push_back(',');
            } else {
                first = false;
            }
            parentIDs.append(clusterID);
        }

        parentClusters.append(std::format("{}:{};", sampleID, parentIDs));
    }

    lineTokens.emplace_back(std::move(parentClusters));

    const std::string color = helper::generateRandomRGBString();
    lineTokens.emplace_back(color);

    writer << lineTokens;
}
}  // namespace

void writeInteractions(const PostprocessData& data,
                       const std::unordered_map<int, std::string>& referenceIndexToIDMap,
                       const InteractionGenerator::Result& clusteringResult) {
    std::ofstream outSuperInteractionsGCT(data.superInteractionsCountsGCTOutPath);
    std::ofstream outSuperInteractionBEDPE(data.superInteractionsBEDPEOutPath);

    if (!outSuperInteractionsGCT.is_open()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not open file: ",
                                                            outSuperInteractionsGCT);
    }

    if (!outSuperInteractionsGCT.is_open()) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Could not open file: ",
                                                            outSuperInteractionBEDPE);
    }

    auto gcsWriter = csv::make_tsv_writer(outSuperInteractionsGCT);
    auto bedpeWriter = csv::make_tsv_writer(outSuperInteractionBEDPE);
    csv::set_decimal_places(2);

    const auto sampleIDs = getSampleIDs(data.samples);

    writeSuperInteractionsGCTHeader(sampleIDs, clusteringResult.superInteractions.size(),
                                    outSuperInteractionsGCT, gcsWriter);
    writeSuperInteractionsBEDPEHeader(outSuperInteractionBEDPE, sampleIDs);

    size_t index = 0;

    for (const auto& superInteraction : clusteringResult.superInteractions) {
        const std::string interactionID = std::format("super_interaction{}", index);

        writeSuperInteractionGCT(gcsWriter, interactionID, sampleIDs, superInteraction);
        writeSuperInteractionsBEDPE(bedpeWriter, referenceIndexToIDMap, interactionID, sampleIDs,
                                    superInteraction);

        ++index;
    }
}

}  // namespace pipelines::postprocess::SuperInteractionsWriter
