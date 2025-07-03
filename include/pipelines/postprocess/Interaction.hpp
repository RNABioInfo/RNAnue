#pragma once

// Standard
#include <algorithm>
#include <cstddef>
#include <ostream>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Internal
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "SortedGenomicRegionPair.hpp"

namespace pipelines::postprocess {
using namespace dataTypes;

struct InteractionID {
    std::string sampleID;
    std::string clusterID;
};

struct InteractionMetrics {
    float contributionScore;
    float meanInterCrosslinkCount;
    float sdInterCrosslinkCount;
    float globalComplementarityScore;
    float globalHybridizationScore;
    float pValue;
    float padjValue;
};

struct InteractionFeatureIDs {
    std::string firstFeature;
    std::string secondFeature;
};

struct InteractionTypeCounts {
    size_t intraInteractionCount;
    size_t interInteractionCount;
};

struct Interaction {
    Interaction(InteractionID interactionID, SortedGenomicRegionPair sortedSegments,
                const InteractionFeatureIDs& interactionFeatureIDs,
                InteractionMetrics interactionMetrics)
        : interactionIDs({std::move(interactionID)}),
          sortedSegments(sortedSegments),
          interactionMetrics({interactionMetrics}),
          firstFeatureIDs({interactionFeatureIDs.firstFeature}),
          secondFeatureIDs({interactionFeatureIDs.secondFeature}),
          transcriptContribution(interactionMetrics.contributionScore) {}

    [[nodiscard]] auto getInteractionIDs() const noexcept -> const std::vector<InteractionID>& {
        return interactionIDs;
    }
    [[nodiscard]] auto getInteractionIDs(std::string_view sampleID) const noexcept {
        return interactionIDs | std::views::filter([&sampleID](const auto& interaction) {
                   return interaction.sampleID == sampleID;
               }) |
               std::views::transform(
                   [](const auto& interaction) -> std::string { return interaction.clusterID; }) |
               std::ranges::to<std::vector>();
    }

    [[nodiscard]] auto getFirstSegment() const noexcept -> const GenomicRegion& {
        return sortedSegments.firstRegion;
    }
    [[nodiscard]] auto getSecondSegment() const noexcept -> const GenomicRegion& {
        return sortedSegments.secondRegion;
    }
    [[nodiscard]] auto getInteractionMetrics() const noexcept
        -> const std::vector<InteractionMetrics>& {
        return interactionMetrics;
    }

    [[nodiscard]] auto getInteractionTypeCounts() const noexcept -> InteractionTypeCounts {
        size_t intraCount{0};
        size_t interCount{0};

        for (const auto& [firstID, secondID] :
             std::ranges::zip_view(firstFeatureIDs, secondFeatureIDs)) {
            if (firstID == secondID) {
                ++intraCount;
            } else {
                ++interCount;
            }
        }

        return {.intraInteractionCount = intraCount, .interInteractionCount = interCount};
    }

    [[nodiscard]] auto getContributionScore(const std::string_view sampleID) const noexcept
        -> float {
        float score = 0;

        for (size_t index : getIndices(sampleID)) {
            score += interactionMetrics[index].contributionScore;
        }

        return score;
    }

    [[nodiscard]] auto getFirstFeatureIDs(const std::string_view sampleID) const noexcept
        -> std::vector<std::string> {
        std::vector<std::string> featureIDs;

        for (size_t index : getIndices(sampleID)) {
            featureIDs.emplace_back(firstFeatureIDs[index]);
        }

        return featureIDs;
    }

    [[nodiscard]] auto getSecondFeatureIDs(const std::string_view sampleID) const noexcept
        -> std::vector<std::string> {
        std::vector<std::string> featureIDs;

        for (size_t index : getIndices(sampleID)) {
            featureIDs.emplace_back(secondFeatureIDs[index]);
        }

        return featureIDs;
    }

    [[nodiscard]] auto getFeatureIDs(const std::string_view sampleID) const noexcept
        -> std::vector<InteractionFeatureIDs> {
        std::vector<InteractionFeatureIDs> featureIDs;

        for (size_t index : getIndices(sampleID)) {
            featureIDs.emplace_back(firstFeatureIDs[index], secondFeatureIDs[index]);
        }

        return featureIDs;
    }

    [[nodiscard]] auto getIndices(const std::string_view sampleID) const noexcept
        -> std::vector<size_t> {
        std::vector<size_t> indices;

        for (size_t index = 0; index < interactionIDs.size(); ++index) {
            if (interactionIDs[index].sampleID == sampleID) {
                indices.emplace_back(index);
            }
        }

        return indices;
    }

    friend auto operator<<(std::ostream& ostream, const Interaction& interaction) -> std::ostream&;

    auto operator<(const Interaction& other) const noexcept -> bool {
        if (getSecondSegment().getReferenceIDIndex() <
            other.getSecondSegment().getReferenceIDIndex()) {
            return true;
        }

        if (getSecondSegment().getReferenceIDIndex() ==
            other.getSecondSegment().getReferenceIDIndex()) {
            return getSecondSegment().getEnd() < other.getSecondSegment().getEnd();
        }
        return false;
    }

    auto operator>(const Interaction& other) const noexcept -> bool { return other < *this; }

    [[nodiscard]] auto isBefore(const Interaction& other) const noexcept -> bool {
        return (getSecondSegment().getReferenceIDIndex() <
                other.getSecondSegment().getReferenceIDIndex()) ||
               (getSecondSegment().getReferenceIDIndex() ==
                    other.getSecondSegment().getReferenceIDIndex() &&
                getSecondSegment().getEnd() < other.getSecondSegment().getStart());
    };

    [[nodiscard]] auto overlapsWithShortestSegmentFraction(
        const Interaction& other, const GenomicStrandSpecificity strandSpecificity,
        float shortestOverlapFraction) const noexcept -> bool {
        return sortedSegments.firstRegion.overlapsWithShortestSegmentFraction(
                   other.getFirstSegment(),
                   GenomicOrientation::fromStrandSpecificity(strandSpecificity),
                   shortestOverlapFraction) &&
               sortedSegments.secondRegion.overlapsWithShortestSegmentFraction(
                   other.getSecondSegment(),
                   GenomicOrientation::fromStrandSpecificity(strandSpecificity),
                   shortestOverlapFraction);
    }

    auto merge(const Interaction& other, const GenomicStrandSpecificity& strandSpecificity) noexcept
        -> bool {
        if (!sortedSegments.merge(other.sortedSegments, strandSpecificity)) [[unlikely]] {
            Logger::log<LogLevel::WARNING>("Could not merge interaction clusters:\n", *this, other);
            return false;
        }

        {
            const auto oldSize = interactionIDs.size();
            const auto otherSize = other.interactionIDs.size();
            interactionIDs.resize(oldSize + otherSize);
            std::ranges::move(other.interactionIDs,
                              interactionIDs.begin() + static_cast<long>(oldSize));
        }

        {
            const auto oldSize = interactionMetrics.size();
            const auto otherSize = other.interactionMetrics.size();
            interactionMetrics.resize(oldSize + otherSize);
            std::ranges::move(other.interactionMetrics,
                              interactionMetrics.begin() + static_cast<long>(oldSize));
        }

        {
            const auto oldSize = firstFeatureIDs.size();
            const auto otherSize = other.firstFeatureIDs.size();
            firstFeatureIDs.resize(oldSize + otherSize);
            std::ranges::move(other.firstFeatureIDs,
                              firstFeatureIDs.begin() + static_cast<long>(oldSize));
        }

        {
            const auto oldSize = secondFeatureIDs.size();
            const auto otherSize = other.secondFeatureIDs.size();
            secondFeatureIDs.resize(oldSize + otherSize);
            std::ranges::move(other.secondFeatureIDs,
                              secondFeatureIDs.begin() + static_cast<long>(oldSize));
        }

        {
            transcriptContribution += other.transcriptContribution;
        }

        return true;
    }

    friend auto operator<<(std::ostream& ostream, const Interaction& interaction) -> std::ostream&;

   private:
    std::vector<InteractionID> interactionIDs;
    SortedGenomicRegionPair sortedSegments;
    std::vector<InteractionMetrics> interactionMetrics;
    std::vector<std::string> firstFeatureIDs;
    std::vector<std::string> secondFeatureIDs;
    float transcriptContribution;
};

inline auto operator<<(std::ostream& ostream, const Interaction& interaction) -> std::ostream& {
    ostream << "Interaction:" << "\n";

    // Print all Interaction IDs
    ostream << "  IDs:" << "\n";
    if (!interaction.interactionIDs.empty()) {
        for (const auto& interactionID : interaction.interactionIDs) {
            ostream << "    { sampleID: " << interactionID.sampleID
                    << ", clusterID: " << interactionID.clusterID << " }"
                    << "\n";
        }
    } else {
        ostream << "    [No IDs available]" << "\n";
    }

    // Print genomic segments using getFirstSegment and getSecondSegment
    ostream << "  Segments:" << "\n";
    ostream << "    First: " << interaction.getFirstSegment() << "\n";
    ostream << "    Second: " << interaction.getSecondSegment() << "\n";

    // Print interaction metrics
    ostream << "  Interaction Metrics:" << "\n";
    if (!interaction.interactionMetrics.empty()) {
        for (const auto& metric : interaction.interactionMetrics) {
            ostream << "    { contributionScore: " << metric.contributionScore
                    << ", meanInterCrosslinkCount: " << metric.meanInterCrosslinkCount
                    << ", sdInterCrosslinkCount: " << metric.sdInterCrosslinkCount
                    << ", globalComplementarityScore: " << metric.globalComplementarityScore
                    << ", globalHybridizationScore: " << metric.globalHybridizationScore
                    << ", pValue: " << metric.pValue << ", padjValue: " << metric.padjValue << " }"
                    << "\n";
        }
    } else {
        ostream << "    [No Metrics available]" << "\n";
    }

    // Print feature IDs
    ostream << "  Feature IDs:" << "\n";
    ostream << "    First Features: ";
    if (!interaction.firstFeatureIDs.empty()) {
        for (const auto& feat : interaction.firstFeatureIDs) {
            ostream << feat << " ";
        }
    } else {
        ostream << "N/A";
    }
    ostream << "\n";

    ostream << "    Second Features: ";
    if (!interaction.secondFeatureIDs.empty()) {
        for (const auto& feat : interaction.secondFeatureIDs) {
            ostream << feat << " ";
        }
    } else {
        ostream << "N/A";
    }
    ostream << "\n";

    // Print transcript contribution
    ostream << "  Transcript Contribution: " << interaction.transcriptContribution << "\n";

    return ostream;
}

}  // namespace pipelines::postprocess
