#pragma once

// Standard
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

// Internal
#include "InteractionSegment.hpp"
#include "RecordFragment.hpp"

namespace pipelines::analyze {

class InteractionCluster {
   public:
    InteractionCluster(InteractionSegmentPair interactionSegments,
                       std::vector<std::string> recordIDs,
                       std::vector<double> complementarityScores,
                       std::vector<double> hybridizationEnergies,
                       std::vector<int32_t> crosslinkingSiteCounts)
        : firstSegment(interactionSegments.firstSegment),
          secondSegment(interactionSegments.secondSegment),
          recordIDs(std::move(recordIDs)),
          complementarityScores(std::move(complementarityScores)),
          maxComplementarityScore(*std::ranges::max_element(this->complementarityScores)),
          hybridizationEnergies(std::move(hybridizationEnergies)),
          minHybridizationEnergy(*std::ranges::min_element(this->hybridizationEnergies)),
          crosslinkingSiteCounts(std::move(crosslinkingSiteCounts)) {}

    InteractionCluster(InteractionSegment firstSegment, InteractionSegment secondSegment,
                       std::string recordID, double complementarityScore,
                       double hybridizationEnergie, int32_t crosslinkingSiteCount)
        : InteractionCluster({.firstSegment = firstSegment, .secondSegment = secondSegment},
                             {std::move(recordID)}, {complementarityScore}, {hybridizationEnergie},
                             {crosslinkingSiteCount}) {}

    InteractionCluster() = delete;

    static auto fromRecordFragments(const RecordFragment &firstFragment,
                                    const RecordFragment &secondFragment) -> InteractionCluster;

    // Getters
    [[nodiscard]] auto getFirstSegment() const -> const InteractionSegment & {
        return firstSegment;
    }

    [[nodiscard]] auto getSecondSegment() const -> const InteractionSegment & {
        return secondSegment;
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
    auto operator<(const InteractionCluster &other) const -> bool;
    auto operator>(const InteractionCluster &other) const -> bool;
    auto operator==(const InteractionCluster &other) const -> bool;

    [[nodiscard]] auto isBefore(const InteractionCluster &other) const noexcept -> bool;

    [[nodiscard]] auto overlaps(const InteractionCluster &other, int graceDistance) const noexcept
        -> bool;

    // Returns the fraction of the overlap between the two segments relative to the length
    // of the shorter segment
    [[nodiscard]] auto segmentsMaxOverlapFraction() const noexcept -> double;

    void merge(const InteractionCluster &other);

    [[nodiscard]] auto complementarityStatistics() const -> double;

    [[nodiscard]] auto hybridizationEnergyStatistics() const -> double;

   private:
    InteractionSegment firstSegment;
    InteractionSegment secondSegment;
    std::vector<std::string> recordIDs;
    std::vector<double> complementarityScores;
    double maxComplementarityScore;
    std::vector<double> hybridizationEnergies;
    double minHybridizationEnergy;
    std::vector<int32_t> crosslinkingSiteCounts;
};

auto operator<<(std::ostream &outputStream, const InteractionCluster &interactionCluster)
    -> std::ostream &;

}  // namespace pipelines::analyze
