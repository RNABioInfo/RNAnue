#pragma once

// Standard
#include <cstddef>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Internal
#include "Constants.hpp"
#include "FeatureParser.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "Interaction.hpp"
#include "Logger.hpp"
#include "Region.hpp"
#include "SortedGenomicRegionPair.hpp"
#include "csv.hpp"

namespace pipelines::postprocess {

namespace fs = std::filesystem;

class InteractionParser {
   public:
    struct Result {
        std::vector<Interaction> interactions;
        std::unordered_map<int, std::string> referenceIndexToIDMap;
    };

    [[nodiscard]] auto getResultAndReset() noexcept -> Result {
        currentIndex = 0;
        referenceIDToIndexMap.clear();
        return {.interactions = std::exchange(interactions, {}),
                .referenceIndexToIDMap = std::exchange(referenceIndexToIDMap, {})};
    }

    void parse(const fs::path& interactionsInputPath, const std::string& sampleID) {
        Logger::log("Parsing interactions file: ", interactionsInputPath);
        csv::CSVFormat format{};
        format.delimiter('\t').header_row(0);
        csv::CSVReader reader{interactionsInputPath.string()};

        using namespace constants::interaction;
        for (const auto& entry : reader) {
            int firstReferenceIndex{
                getReferenceIndexForID(entry[firstSegmentReferenceHeader].get<std::string>())};
            GenomicRegion firstRegion{
                firstReferenceIndex,
                Region{.startPosition = entry[firstSegmentStartHeader].get<int>(),
                       .endPosition = entry[firstSegmentEndHeader].get<int>()},
                getGenomicStrand(*entry[firstSegmentStrandHeader].get<std::string>().cbegin())};

            int secondReferenceIndex{
                getReferenceIndexForID(entry[secondSegmentReferenceHeader].get<std::string>())};
            GenomicRegion secondRegion{
                secondReferenceIndex,
                Region{.startPosition = entry[secondSegmentStartHeader].get<int>(),
                       .endPosition = entry[secondSegmentEndHeader].get<int>()},
                getGenomicStrand(*entry[secondSegmentStrandHeader].get<std::string>().cbegin())};

            interactions.emplace_back(
                InteractionID{.sampleID = sampleID,
                              .clusterID = entry[clusterIDHeader].get<std::string>()},
                SortedGenomicRegionPair(firstRegion, secondRegion),
                InteractionFeatureIDs{
                    .firstFeature = entry[firstFeatureIDHeader].get<std::string>(),
                    .secondFeature = entry[secondFeatureIDHeader].get<std::string>()},
                InteractionMetrics{
                    .contributionScore = entry[transcriptContributionHeader].get<float>(),
                    .meanInterCrosslinkCount = entry[meanInterCrosslinkCountHeader].get<float>(),
                    .sdInterCrosslinkCount = entry[sdInterCrosslinksHeader].get<float>(),
                    .globalComplementarityScore =
                        entry[globalHybridizationScoreHeader].get<float>(),
                    .globalHybridizationScore =
                        entry[globalComplementarityScoreHeader].get<float>(),
                    .pValue = entry[pValueHeader].get<float>(),
                    .padjValue = entry[padjValueHeader].get<float>()});
        }
    }

   private:
    annotation::ReferenceIDToIndexMap referenceIDToIndexMap;
    std::unordered_map<int, std::string> referenceIndexToIDMap;
    int currentIndex{0};
    std::vector<Interaction> interactions;

    [[nodiscard]] auto getReferenceIndexForID(const std::string& referenceID) -> int {
        if (!referenceIDToIndexMap.contains(referenceID)) {
            int lastIndex = currentIndex;

            referenceIDToIndexMap.emplace(referenceID, lastIndex);
            referenceIndexToIDMap.emplace(lastIndex, referenceID);
            ++currentIndex;

            return lastIndex;
        }

        return referenceIDToIndexMap[referenceID];
    }
};

}  // namespace pipelines::postprocess
