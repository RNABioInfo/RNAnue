#pragma once

// Standard
#include <algorithm>
#include <cstddef>
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
    InteractionCluster(InteractionSegment firstSegment, InteractionSegment secondSegment,
                       std::vector<std::string> firstSegmentRecordIDs,
                       std::vector<std::string> secondSegmentRecordIDs,
                       std::vector<double> complementarityScores,
                       std::vector<double> hybridizationEnergies)
        : firstSegment(firstSegment),
          secondSegment(secondSegment),
          firstSegmentRecordIDs(std::move(firstSegmentRecordIDs)),
          secondSegmentRecordIDs(std::move(secondSegmentRecordIDs)),
          complementarityScores(std::move(complementarityScores)),
          maxComplementarityScore(*std::ranges::max_element(this->complementarityScores)),
          hybridizationEnergies(std::move(hybridizationEnergies)),
          minHybridizationEnergy(*std::ranges::min_element(this->hybridizationEnergies)) {}

    InteractionCluster(InteractionSegment firstSegment, InteractionSegment secondSegment,
                       std::string firstSegmentRecordID, std::string secondSegmentRecordID,
                       double complementarityScore, double hybridizationEnergie)
        : InteractionCluster(firstSegment, secondSegment,
                             std::vector<std::string>{std::move(firstSegmentRecordID)},
                             std::vector<std::string>{std::move(secondSegmentRecordID)},
                             std::vector<double>{complementarityScore},
                             std::vector<double>{hybridizationEnergie}) {}

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

    [[nodiscard]] auto getFirstSegmentRecordIDs() const -> const std::vector<std::string> & {
        return firstSegmentRecordIDs;
    }

    [[nodiscard]] auto getSecondSegmentRecordIDs() const -> const std::vector<std::string> & {
        return secondSegmentRecordIDs;
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

    [[nodiscard]] auto fragmentCount() const -> size_t { return firstSegmentRecordIDs.size(); }

    auto operator<(const InteractionCluster &other) const -> bool;
    auto operator>(const InteractionCluster &other) const -> bool;
    auto operator==(const InteractionCluster &other) const -> bool;

    [[nodiscard]] auto overlaps(const InteractionCluster &other, int graceDistance) const noexcept
        -> bool;

    void merge(const InteractionCluster &other);

    [[nodiscard]] auto complementarityStatistics() const -> double;

    [[nodiscard]] auto hybridizationEnergyStatistics() const -> double;

   private:
    InteractionSegment firstSegment;
    InteractionSegment secondSegment;
    std::vector<std::string> firstSegmentRecordIDs;
    std::vector<std::string> secondSegmentRecordIDs;
    std::vector<double> complementarityScores;
    double maxComplementarityScore;
    std::vector<double> hybridizationEnergies;
    double minHybridizationEnergy;
};

auto operator<<(std::ostream &outputStream, const InteractionCluster &interactionCluster)
    -> std::ostream &;

}  // namespace pipelines::analyze
