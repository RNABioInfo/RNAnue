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

    return {firstSegment,
            secondSegment,
            firstFragment.recordID,
            secondFragment.recordID,
            firstFragment.complementarityScore,
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
    } else if (secondSegment.getReferenceIDIndex() == other.secondSegment.getReferenceIDIndex()) {
        return secondSegment.getEnd() < other.secondSegment.getEnd();
    } else {
        return false;
    }
}

auto InteractionCluster::operator>(const InteractionCluster &other) const -> bool {
    return other < *this;
}

auto InteractionCluster::operator==(const InteractionCluster &other) const -> bool {
    return firstSegment == other.firstSegment && secondSegment == other.secondSegment &&
           firstSegmentRecordIDs == other.firstSegmentRecordIDs &&
           secondSegmentRecordIDs == other.secondSegmentRecordIDs &&
           complementarityScores == other.complementarityScores &&
           hybridizationEnergies == other.hybridizationEnergies;
}

auto InteractionCluster::overlaps(const InteractionCluster &other,
                                  const int graceDistance) const noexcept -> bool {
    return firstSegment.overlaps(other.firstSegment, graceDistance) &&
           secondSegment.overlaps(other.secondSegment, graceDistance);
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

    firstSegmentRecordIDs.reserve(firstSegmentRecordIDs.size() +
                                  other.firstSegmentRecordIDs.size());
    firstSegmentRecordIDs.insert(firstSegmentRecordIDs.end(),
                                 std::make_move_iterator(other.firstSegmentRecordIDs.begin()),
                                 std::make_move_iterator(other.firstSegmentRecordIDs.end()));

    secondSegmentRecordIDs.reserve(secondSegmentRecordIDs.size() +
                                   other.secondSegmentRecordIDs.size());
    secondSegmentRecordIDs.insert(secondSegmentRecordIDs.end(),
                                  std::make_move_iterator(other.secondSegmentRecordIDs.begin()),
                                  std::make_move_iterator(other.secondSegmentRecordIDs.end()));

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
    return outputStream << "InteractionCluster:\n"
                        << "First segment: " << interactionCluster.getFirstSegment()
                        << "\nSecond segment id: " << interactionCluster.getSecondSegment() << "\n"
                        << "Complementarity score count: "
                        << interactionCluster.getComplementarityScores().size() << "\n"
                        << "Hybridization energie count: "
                        << interactionCluster.getHybridizationEnergies().size() << "\n"
                        << "Count: " << interactionCluster.fragmentCount() << "\n";
}

}  // namespace pipelines::analyze
