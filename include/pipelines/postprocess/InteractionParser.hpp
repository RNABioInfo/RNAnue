#pragma once

// Standard
#include <cstddef>
#include <filesystem>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
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
        const auto colNames = reader.get_col_names();
        const std::unordered_set<std::string> columns{colNames.begin(), colNames.end()};

        using namespace constants::interaction;
        const auto hasColumn = [&columns](const std::string& columnName) -> bool {
            return columns.contains(columnName);
        };
        const auto optionalDouble = [&hasColumn](const csv::CSVRow& entry,
                                                const std::string& columnName) -> double {
            return hasColumn(columnName) ? entry[columnName].get<double>()
                                         : std::numeric_limits<double>::quiet_NaN();
        };
        const auto optionalString = [&hasColumn](const csv::CSVRow& entry,
                                                 const std::string& columnName) -> std::string {
            return hasColumn(columnName) ? entry[columnName].get<std::string>() : std::string{};
        };

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
                    .contributionScore = entry[transcriptContributionHeader].get<double>(),
                    .totalSpanBp = optionalDouble(entry, totalSpanBpHeader),
                    .effectiveCoverageSpanBp =
                        optionalDouble(entry, effectiveCoverageSpanBpHeader),
                    .supportPerTotalBp = optionalDouble(entry, supportPerTotalBpHeader),
                    .supportPerEffectiveBp = optionalDouble(entry, supportPerEffectiveBpHeader),
                    .coverageConcentration = optionalDouble(entry, coverageConcentrationHeader),
                    .coverageComponents = optionalDouble(entry, coverageComponentsHeader),
                    .armBalance = optionalDouble(entry, armBalanceHeader),
                    .coverageProfile = optionalString(entry, coverageProfileHeader),
                    .meanInterCrosslinkCount = entry[meanInterCrosslinkCountHeader].get<double>(),
                    .sdInterCrosslinkCount = entry[sdInterCrosslinksHeader].get<double>(),
                    .globalComplementarityScore =
                        entry[globalComplementarityScoreHeader].get<double>(),
                    .globalHybridizationScore =
                        entry[globalHybridizationScoreHeader].get<double>(),
                    .pValue = entry[pValueHeader].get<double>(),
                    .padjValue = entry[padjValueHeader].get<double>()});
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
