#pragma once

// Standard
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

// Internal
#include "ArmCoverage.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "RecordFragment.hpp"
#include "SortedGenomicRegionPair.hpp"

namespace pipelines::analyze {

using namespace dataTypes;

struct CoverageShapeMetrics {
    size_t totalSpanBp{};
    double effectiveCoverageSpanBp{};
    double supportPerTotalBp{};
    double supportPerEffectiveBp{};
    double coverageConcentration{};
    size_t coverageComponents{};
    double armBalance{};
    std::string coverageProfile;
};

class InteractionCluster {
   public:
    InteractionCluster(SortedGenomicRegionPair sortedSegments, std::vector<std::string> recordIDs,
                       std::vector<double> complementarityScores,
                       std::vector<double> hybridizationEnergies,
                       std::vector<int32_t> crosslinkingSiteCounts, float transcriptContribution)
        : sortedSegments(sortedSegments),
          recordIDs(std::move(recordIDs)),
          complementarityScores(std::move(complementarityScores)),
          maxComplementarityScore(*std::ranges::max_element(this->complementarityScores)),
          hybridizationEnergies(std::move(hybridizationEnergies)),
          minHybridizationEnergy(*std::ranges::min_element(this->hybridizationEnergies)),
          crosslinkingSiteCounts(std::move(crosslinkingSiteCounts)),
          transcriptContribution(transcriptContribution),
          transcriptContributionSquaredSum(estimateContributionSquareSum(
              this->recordIDs.size(), transcriptContribution)) {
        addFullSpanCoverage(transcriptContribution);
    };

    InteractionCluster(SortedGenomicRegionPair sortedSegments, std::string recordID,
                       double complementarityScore, double hybridizationEnergy,
                       int crosslinkingSiteCount, float transcriptContribution)
        : sortedSegments(sortedSegments),
          recordIDs({std::move(recordID)}),
          complementarityScores({complementarityScore}),
          maxComplementarityScore(complementarityScore),
          hybridizationEnergies({hybridizationEnergy}),
          minHybridizationEnergy(hybridizationEnergy),
          crosslinkingSiteCounts({crosslinkingSiteCount}),
          transcriptContribution(transcriptContribution),
          transcriptContributionSquaredSum(transcriptContribution * transcriptContribution) {
        addFullSpanCoverage(transcriptContribution);
    };

    InteractionCluster() = delete;

    static auto fromRecordFragments(const RecordFragment &firstFragment,
                                    const RecordFragment &secondFragment)
        -> InteractionCluster {
        assert((firstFragment.recordID == secondFragment.recordID) &&
               (firstFragment.complementarityScore == secondFragment.complementarityScore) &&
               (firstFragment.hybridizationEnergy == secondFragment.hybridizationEnergy) &&
               (firstFragment.interCrosslinkingSiteCount ==
                secondFragment.interCrosslinkingSiteCount) &&
               (firstFragment.transcriptContribution == secondFragment.transcriptContribution));

        InteractionCluster cluster{{firstFragment.genomicRegion, secondFragment.genomicRegion},
                                   firstFragment.recordID,
                                   firstFragment.complementarityScore,
                                   firstFragment.hybridizationEnergy,
                                   firstFragment.interCrosslinkingSiteCount,
                                   firstFragment.transcriptContribution};
        cluster.replaceCoverageFromFragments(firstFragment, secondFragment);
        return cluster;
    };

    // Getters
    [[nodiscard]] auto getFirstSegment() const -> const GenomicRegion & {
        return sortedSegments.firstRegion;
    }

    [[nodiscard]] auto getSecondSegment() const -> const GenomicRegion & {
        return sortedSegments.secondRegion;
    }

    [[nodiscard]] auto getRecordIDs() const -> const std::vector<std::string> & {
        return recordIDs;
    }

    [[nodiscard]] auto fragmentCount() const noexcept -> size_t { return recordIDs.size(); }

    [[nodiscard]] auto getComplementarityScores() const -> const std::vector<double> & {
        return complementarityScores;
    }

    [[nodiscard]] auto getMaxComplementarityScore() const -> double {
        return maxComplementarityScore;
    }

    [[nodiscard]] auto getHybridizationEnergies() const -> const std::vector<double> & {
        return hybridizationEnergies;
    }

    [[nodiscard]] auto getMinHybridizationEnergy() const -> double {
        return minHybridizationEnergy;
    }

    [[nodiscard]] auto getCrosslinkingSiteCounts() const -> const std::vector<int32_t> & {
        return crosslinkingSiteCounts;
    }

    [[nodiscard]] auto meanCrosslinkingSiteCount() const -> double;

    [[nodiscard]] auto standardDeviationCrosslinkingSiteCount() const -> double;

    [[nodiscard]] auto getTranscriptContribution() const -> float { return transcriptContribution; }

    [[nodiscard]] auto getTranscriptContributionSquaredSum() const -> double {
        return transcriptContributionSquaredSum;
    }

    [[nodiscard]] auto totalSpanBp() const noexcept -> size_t {
        return getFirstSegment().length() + getSecondSegment().length();
    }

    [[nodiscard]] auto maxArmSpanBp() const noexcept -> size_t {
        return std::max(getFirstSegment().length(), getSecondSegment().length());
    }

    [[nodiscard]] auto coverageShapeMetrics() const -> CoverageShapeMetrics;

    [[nodiscard]] auto getFirstArmCoverageRuns() const -> std::vector<CoverageRun> {
        return firstArmCoverage.runs();
    }

    [[nodiscard]] auto getSecondArmCoverageRuns() const -> std::vector<CoverageRun> {
        return secondArmCoverage.runs();
    }

    // Comparisons
    /**
     * @brief Less-than comparison operator for InteractionCluster.
     *
     * This operator compares two InteractionCluster objects based on the end position
     * of the second segment and the referenceIndexID. It returns true if the referenceIndexID is
     * less in the current object or the end position of the segment in the current object are
     * lexicographically less than those in the provided object.
     *
     * @param a The InteractionCluster object to compare with.
     * @return true if the current object is less than the provided object, false otherwise.
     */
    auto operator<(const InteractionCluster &other) const noexcept -> bool;
    auto operator>(const InteractionCluster &other) const noexcept -> bool;
    auto operator==(const InteractionCluster &other) const noexcept -> bool;

    /**
     * @brief Determines if this InteractionCluster is strictly located before another
     * InteractionCluster.
     *
     * This method compares the second segment of the current InteractionCluster with the second
     * segment of the provided one based on their reference index and positions to check if it
     * appears entirely before, with no overlap and no blunt end.
     *
     * @param other The InteractionCluster to compare against.
     * @return true if this InteractionCluster is entirely before the other, false otherwise.
     */
    [[nodiscard]] auto isBefore(const InteractionCluster &other) const noexcept -> bool;

    [[nodiscard]] inline auto overlapsWithTolerance(const InteractionCluster &other,
                                                    GenomicStrandSpecificity strandSpecificity,
                                                    int tolerance) const noexcept -> bool {
        return sortedSegments.firstRegion.overlapsWithTolerance(
                   other.getFirstSegment(), overlapOrientation(strandSpecificity), tolerance) &&
               sortedSegments.secondRegion.overlapsWithTolerance(
                   other.getSecondSegment(), overlapOrientation(strandSpecificity), tolerance);
    };

    [[nodiscard]] inline auto overlapsWithShortestSegmentFraction(
        const InteractionCluster &other, GenomicStrandSpecificity strandSpecificity,
        float shortestOverlapFraction) const noexcept -> bool {
        return sortedSegments.firstRegion.overlapsWithShortestSegmentFraction(
                   other.getFirstSegment(),
                   GenomicOrientation::fromStrandSpecificity(strandSpecificity),
                   shortestOverlapFraction) &&
               sortedSegments.secondRegion.overlapsWithShortestSegmentFraction(
                   other.getSecondSegment(),
                   GenomicOrientation::fromStrandSpecificity(strandSpecificity),
                   shortestOverlapFraction);
    };

    // Returns the fraction of the overlap between the two segments relative to the length
    // of the shorter segment
    [[nodiscard]] auto segmentsMaxSelfOverlapFraction() const noexcept -> double;

    [[nodiscard]] auto merge(const InteractionCluster &other,
                             const GenomicStrandSpecificity &) noexcept -> bool;

    void absorbValidatedComponentMember(const InteractionCluster &other) noexcept;

    [[nodiscard]] auto complementarityStatistics() const -> double;

    [[nodiscard]] auto hybridizationEnergyStatistics() const -> double;

   private:
    SortedGenomicRegionPair sortedSegments;
    std::vector<std::string> recordIDs;
    std::vector<double> complementarityScores;
    double maxComplementarityScore;
    std::vector<double> hybridizationEnergies;
    double minHybridizationEnergy;
    std::vector<int32_t> crosslinkingSiteCounts;
    float transcriptContribution;
    double transcriptContributionSquaredSum;
    ArmCoverage firstArmCoverage;
    ArmCoverage secondArmCoverage;

    [[nodiscard]] static constexpr auto estimateContributionSquareSum(
        size_t recordCount, float transcriptContribution) noexcept -> double {
        if (recordCount == 0) {
            return 0.0;
        }

        const double uniformContribution =
            static_cast<double>(transcriptContribution) / static_cast<double>(recordCount);
        return static_cast<double>(recordCount) * uniformContribution * uniformContribution;
    }

    [[nodiscard]] static constexpr auto overlapOrientation(
        const GenomicStrandSpecificity strandSpecificity) -> GenomicOrientation {
        return (strandSpecificity == GenomicStrandSpecificity::SPECIFIC) ? GenomicOrientation::SAME
                                                                         : GenomicOrientation::BOTH;
    };

    void addFullSpanCoverage(double weight) {
        firstArmCoverage.addInterval(getFirstSegment().getStart(), getFirstSegment().getEnd(),
                                     weight);
        secondArmCoverage.addInterval(getSecondSegment().getStart(), getSecondSegment().getEnd(),
                                      weight);
    }

    void replaceCoverageFromFragments(const RecordFragment &firstFragment,
                                      const RecordFragment &secondFragment);
};

auto operator<<(std::ostream &outputStream, const InteractionCluster &interactionCluster)
    -> std::ostream &;

}  // namespace pipelines::analyze
