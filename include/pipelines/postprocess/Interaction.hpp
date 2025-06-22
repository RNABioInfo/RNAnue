#pragma once

// Standard
#include <algorithm>
#include <ostream>
#include <string>
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
            interactionMetrics.resize(oldSize + otherSize);
            std::ranges::move(other.firstFeatureIDs,
                              firstFeatureIDs.begin() + static_cast<long>(oldSize));
        }

        {
            const auto oldSize = secondFeatureIDs.size();
            const auto otherSize = other.secondFeatureIDs.size();
            interactionMetrics.resize(oldSize + otherSize);
            std::ranges::move(other.secondFeatureIDs,
                              firstFeatureIDs.begin() + static_cast<long>(oldSize));
        }

        {
            transcriptContribution += other.transcriptContribution;
        }

        return true;
    }

   private:
    std::vector<InteractionID> interactionIDs;
    SortedGenomicRegionPair sortedSegments;
    std::vector<InteractionMetrics> interactionMetrics;
    std::vector<std::string> firstFeatureIDs;
    std::vector<std::string> secondFeatureIDs;
    float transcriptContribution;
};

inline auto operator<<(std::ostream& ostream, const Interaction& interaction) -> std::ostream& {
    ostream << "Interaction ";
    if (!interaction.interactionIDs.empty()) {
        ostream << "[sampleID: " << interaction.interactionIDs.front().sampleID
                << ", clusterID: " << interaction.interactionIDs.front().clusterID << "] ";
    } else {
        ostream << "[no ID] ";
    }
    ostream << "Segments { First: " << interaction.getFirstSegment()
            << ", Second: " << interaction.getSecondSegment() << " } ";

    if (!interaction.interactionMetrics.empty()) {
        const auto& metrics = interaction.interactionMetrics.front();
        ostream << "Metrics { meanInterCrosslinkCount: " << metrics.meanInterCrosslinkCount
                << ", contributionScore: " << metrics.contributionScore
                << ", pValue: " << metrics.pValue << " } ";
    }

    if (!interaction.firstFeatureIDs.empty() || !interaction.secondFeatureIDs.empty()) {
        ostream << "Features { First: ";
        ostream << (!interaction.firstFeatureIDs.empty() ? interaction.firstFeatureIDs.front()
                                                         : "N/A");
        ostream << ", Second: ";
        ostream << (!interaction.secondFeatureIDs.empty() ? interaction.secondFeatureIDs.front()
                                                          : "N/A");
        ostream << " } ";
    }

    ostream << "TranscriptContribution: " << interaction.transcriptContribution;
    return ostream;
}

}  // namespace pipelines::postprocess
