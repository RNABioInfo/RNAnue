#include "InteractionCluster.hpp"

#include <algorithm>
#include <ostream>

#include "InteractionSegment.hpp"
#include "Utility.hpp"

namespace pipelines::analyze {

auto InteractionCluster::fromRecordFragments(const RecordFragment &firstFragment,
                                             const RecordFragment &secondFragment)
    -> InteractionCluster {
    assert(firstFragment.recordID == secondFragment.recordID);
    assert(firstFragment.complementarityScore == secondFragment.complementarityScore);
    assert(firstFragment.hybridizationEnergy == secondFragment.hybridizationEnergy);

    InteractionSegment firstSegment = {std::min(firstFragment, secondFragment)};
    InteractionSegment secondSegment = {std::max(firstFragment, secondFragment)};

    return {firstSegment, secondSegment, firstFragment.recordID, firstFragment.complementarityScore,
            firstFragment.hybridizationEnergy};
}

/**
 * @brief Less-than comparison operator for InteractionCluster.
 *
 * This operator compares two InteractionCluster objects based on the end position
 * of the second segment and the referenceIndexID. It returns true if the referenceIndexID is less
 * in the current object or the end position of the segment in the current object are
 * lexicographically less than those in the provided object.
 *
 * @param a The InteractionCluster object to compare with.
 * @return true if the current object is less than the provided object, false otherwise.
 */
auto InteractionCluster::operator<(const InteractionCluster &other) const -> bool {
    if (secondSegment.getReferenceIDIndex() < other.secondSegment.getReferenceIDIndex()) {
        return true;
    }
    if (secondSegment.getReferenceIDIndex() == other.secondSegment.getReferenceIDIndex()) {
        return secondSegment.getEnd() < other.secondSegment.getEnd();
    }
    return false;
}

auto InteractionCluster::operator>(const InteractionCluster &other) const -> bool {
    return other < *this;
}

auto InteractionCluster::operator==(const InteractionCluster &other) const -> bool {
    return firstSegment == other.firstSegment && secondSegment == other.secondSegment &&
           recordIDs == other.recordIDs &&
           helper::vectorsApproxEqual(complementarityScores, other.complementarityScores) &&
           helper::vectorsApproxEqual(hybridizationEnergies, other.hybridizationEnergies);
}

auto InteractionCluster::isBefore(const InteractionCluster &other) const noexcept -> bool {
    return (getSecondSegment().getReferenceIDIndex() <
            other.getSecondSegment().getReferenceIDIndex()) ||
           (getSecondSegment().getEnd() < other.getSecondSegment().getStart());
};

auto InteractionCluster::overlaps(const InteractionCluster &other,
                                  const int graceDistance) const noexcept -> bool {
    return firstSegment.overlaps(other.firstSegment, graceDistance) &&
           secondSegment.overlaps(other.secondSegment, graceDistance);
}

auto InteractionCluster::segmentsMaxOverlapFraction() const noexcept -> double {
    return std::max(firstSegment.overlapFraction(secondSegment),
                    secondSegment.overlapFraction(firstSegment));
}

void InteractionCluster::merge(const InteractionCluster &other) {
    assert(firstSegment.getReferenceIDIndex() == other.firstSegment.getReferenceIDIndex());
    assert(secondSegment.getReferenceIDIndex() == other.secondSegment.getReferenceIDIndex());
    assert(firstSegment.getStrand() == other.firstSegment.getStrand());
    assert(secondSegment.getStrand() == other.secondSegment.getStrand());

    firstSegment.merge(other.firstSegment);
    secondSegment.merge(other.secondSegment);

    maxComplementarityScore = std::max(maxComplementarityScore, other.maxComplementarityScore);
    minHybridizationEnergy = std::min(minHybridizationEnergy, other.minHybridizationEnergy);

    recordIDs.reserve(recordIDs.size() + other.recordIDs.size());
    recordIDs.insert(recordIDs.end(), std::make_move_iterator(other.recordIDs.begin()),
                     std::make_move_iterator(other.recordIDs.end()));

    complementarityScores.reserve(complementarityScores.size() +
                                  other.complementarityScores.size());
    complementarityScores.insert(complementarityScores.end(),
                                 std::make_move_iterator(other.complementarityScores.begin()),
                                 std::make_move_iterator(other.complementarityScores.end()));

    hybridizationEnergies.reserve(hybridizationEnergies.size() +
                                  other.hybridizationEnergies.size());
    hybridizationEnergies.insert(hybridizationEnergies.end(),
                                 std::make_move_iterator(other.hybridizationEnergies.begin()),
                                 std::make_move_iterator(other.hybridizationEnergies.end()));
}

auto InteractionCluster::complementarityStatistics() const -> double {
    return helper::calculateMedian(complementarityScores) * maxComplementarityScore;
}

auto InteractionCluster::hybridizationEnergyStatistics() const -> double {
    return std::sqrt(helper::calculateMedian(hybridizationEnergies) * minHybridizationEnergy);
}

auto operator<<(std::ostream &outputStream, const InteractionCluster &interactionCluster)
    -> std::ostream & {
    outputStream << "InteractionCluster:\n"
                 << "First segment: " << interactionCluster.getFirstSegment()
                 << "\nSecond segment id: " << interactionCluster.getSecondSegment() << "\n"
                 << "Complementarity score count: "
                 << interactionCluster.getComplementarityScores().size() << "\n"
                 << "Hybridization energie count: "
                 << interactionCluster.getHybridizationEnergies().size() << "\n"
                 << "Count: " << interactionCluster.fragmentCount() << "\n";

    for (const auto &recordID : interactionCluster.getRecordIDs()) {
        outputStream << recordID << ", ";
    }
    outputStream << "\n";

    return outputStream;
}

}  // namespace pipelines::analyze
