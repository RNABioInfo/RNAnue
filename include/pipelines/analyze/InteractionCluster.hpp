#pragma once

// Standard
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

// Internal
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "RecordFragment.hpp"
#include "SortedGenomicRegionPair.hpp"

namespace pipelines::analyze {

using namespace dataTypes;

class InteractionCluster {
   public:
    InteractionCluster(SortedGenomicRegionPair sortedSegments, std::vector<std::string> recordIDs,
                       std::vector<double> complementarityScores,
                       std::vector<double> hybridizationEnergies,
                       std::vector<int32_t> crosslinkingSiteCounts);

    InteractionCluster(SortedGenomicRegionPair sortedSegments, std::string recordID,
                       double complementarityScore, double hybridizationEnergy,
                       int crosslinkingSiteCount);

    InteractionCluster() = delete;

    static auto fromRecordFragments(const RecordFragment &firstFragment,
                                    const RecordFragment &secondFragment) -> InteractionCluster;

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

    [[nodiscard]] auto fragmentCount() const -> size_t { return recordIDs.size(); }

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

    [[nodiscard]] auto overlapsWithTolerance(const InteractionCluster &other,
                                             GenomicStrandSpecificity strandSpecificity,
                                             int tolerance) const noexcept -> bool;

    [[nodiscard]] auto overlapsWithShortestSegmentFraction(
        const InteractionCluster &other, GenomicStrandSpecificity strandSpecificity,
        float shortestOverlapFraction) const noexcept -> bool;

    // Returns the fraction of the overlap between the two segments relative to the length
    // of the shorter segment
    [[nodiscard]] auto segmentsMaxSelfOverlapFraction() const noexcept -> double;

    [[nodiscard]] auto merge(const InteractionCluster &other,
                             const GenomicStrandSpecificity &) noexcept -> bool;

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

    [[nodiscard]] static constexpr auto overlapOrientation(
        const GenomicStrandSpecificity strandSpecificity) -> GenomicOrientation {
        return (strandSpecificity == GenomicStrandSpecificity::SPECIFIC) ? GenomicOrientation::SAME
                                                                         : GenomicOrientation::BOTH;
    };
};

auto operator<<(std::ostream &outputStream, const InteractionCluster &interactionCluster)
    -> std::ostream &;

}  // namespace pipelines::analyze
