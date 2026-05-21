#include "InteractionCluster.hpp"

// Standard
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ios>
#include <limits>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>

// Internal
#include "GenomicOrientation.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "SortedGenomicRegionPair.hpp"
#include "Utility.hpp"

namespace pipelines::analyze {
namespace {

constexpr double defaultMinSupportPerEffectiveBp = 0.02;
constexpr double minBalancedArmRatio = 0.25;
constexpr double diffuseConcentrationMax = 0.25;
constexpr size_t compactTotalSpanMaxBp = 300;
constexpr size_t expectedSingleComponentCount = 2;

[[nodiscard]] auto safeDivide(double numerator, double denominator) noexcept -> double {
    if (!std::isfinite(numerator) || !std::isfinite(denominator) || denominator <= 0.0) {
        return 0.0;
    }

    return numerator / denominator;
}

[[nodiscard]] auto clampUnitInterval(double value) noexcept -> double {
    if (!std::isfinite(value)) {
        return 0.0;
    }

    return std::clamp(value, 0.0, 1.0);
}

[[nodiscard]] auto fallbackIntervalForRegion(const GenomicRegion &region)
    -> std::vector<CoverageInterval> {
    return {{.start = region.getStart(), .end = region.getEnd()}};
}

[[nodiscard]] auto coverageIntervalsForFragment(const RecordFragment &fragment)
    -> std::vector<CoverageInterval> {
    if (!fragment.coverageIntervals.empty()) {
        return fragment.coverageIntervals;
    }

    return fallbackIntervalForRegion(fragment.genomicRegion);
}

void addFragmentCoverage(ArmCoverage &coverage, const RecordFragment &fragment) {
    coverage.addIntervals(coverageIntervalsForFragment(fragment), fragment.transcriptContribution);
}

[[nodiscard]] auto classifyCoverageProfile(const CoverageShapeMetrics &metrics) -> std::string {
    if (metrics.effectiveCoverageSpanBp <= 0.0 ||
        metrics.supportPerEffectiveBp < defaultMinSupportPerEffectiveBp) {
        return "low_support_density";
    }

    if (metrics.armBalance < minBalancedArmRatio) {
        return "imbalanced_support";
    }

    if (metrics.coverageComponents > expectedSingleComponentCount) {
        return "multi_peak_refine";
    }

    if (metrics.totalSpanBp > compactTotalSpanMaxBp &&
        metrics.coverageConcentration <= diffuseConcentrationMax) {
        return "broad_diffuse";
    }

    return "compact_dense";
}

}  // namespace

auto InteractionCluster::operator<(const InteractionCluster &other) const noexcept -> bool {
    if (getSecondSegment().getReferenceIDIndex() < other.getSecondSegment().getReferenceIDIndex()) {
        return true;
    }

    if (getSecondSegment().getReferenceIDIndex() ==
        other.getSecondSegment().getReferenceIDIndex()) {
        return getSecondSegment().getEnd() < other.getSecondSegment().getEnd();
    }
    return false;
}

auto InteractionCluster::operator>(const InteractionCluster &other) const noexcept -> bool {
    return other < *this;
}

auto InteractionCluster::operator==(const InteractionCluster &other) const noexcept -> bool {
    return sortedSegments == other.sortedSegments && recordIDs == other.recordIDs &&
           helper::vectorsApproxEqual(complementarityScores, other.complementarityScores) &&
           helper::vectorsApproxEqual(hybridizationEnergies, other.hybridizationEnergies) &&
           crosslinkingSiteCounts == other.crosslinkingSiteCounts;
}

auto InteractionCluster::isBefore(const InteractionCluster &other) const noexcept -> bool {
    return (getSecondSegment().getReferenceIDIndex() <
            other.getSecondSegment().getReferenceIDIndex()) ||
           (getSecondSegment().getReferenceIDIndex() ==
                other.getSecondSegment().getReferenceIDIndex() &&
            getSecondSegment().getEnd() < other.getSecondSegment().getStart());
};

auto InteractionCluster::segmentsMaxSelfOverlapFraction() const noexcept -> double {
    return std::max(
        getFirstSegment().overlapFraction(getSecondSegment(), GenomicOrientation::SAME),
        getSecondSegment().overlapFraction(getFirstSegment(), GenomicOrientation::SAME));
}

auto InteractionCluster::merge(const InteractionCluster &other,
                               const GenomicStrandSpecificity &strandSpecificity) noexcept -> bool {
    if (!sortedSegments.merge(other.sortedSegments, strandSpecificity)) [[unlikely]] {
        Logger::log<LogLevel::WARNING>("Could not merge interaction clusters:\n", *this, other);
        return false;
    }

    maxComplementarityScore = (std::max)(maxComplementarityScore, other.maxComplementarityScore);
    minHybridizationEnergy = (std::min)(minHybridizationEnergy, other.minHybridizationEnergy);

    {
        const auto oldSize = recordIDs.size();
        const auto otherSize = other.recordIDs.size();
        recordIDs.resize(oldSize + otherSize);
        std::ranges::move(other.recordIDs, recordIDs.begin() + static_cast<long>(oldSize));
    }

    {
        const auto oldSize = complementarityScores.size();
        const auto otherSize = other.complementarityScores.size();
        complementarityScores.resize(oldSize + otherSize);
        std::ranges::move(other.complementarityScores,
                          complementarityScores.begin() + static_cast<long>(oldSize));
    }

    {
        const auto oldSize = hybridizationEnergies.size();
        const auto otherSize = other.hybridizationEnergies.size();
        hybridizationEnergies.resize(oldSize + otherSize);
        std::ranges::move(other.hybridizationEnergies,
                          hybridizationEnergies.begin() + static_cast<long>(oldSize));
    }

    {
        const auto oldSize = crosslinkingSiteCounts.size();
        const auto otherSize = other.crosslinkingSiteCounts.size();
        crosslinkingSiteCounts.resize(oldSize + otherSize);
        std::ranges::move(other.crosslinkingSiteCounts,
                          crosslinkingSiteCounts.begin() + static_cast<long>(oldSize));
    }

    {
        transcriptContribution += other.transcriptContribution;
        transcriptContributionSquaredSum += other.transcriptContributionSquaredSum;
    }

    firstArmCoverage.merge(other.firstArmCoverage);
    secondArmCoverage.merge(other.secondArmCoverage);

    return true;
}

void InteractionCluster::absorbValidatedComponentMember(const InteractionCluster &other) noexcept {
    auto absorbSegment = [](GenomicRegion &segment, const GenomicRegion &otherSegment) {
        assert(segment.getReferenceIDIndex() == otherSegment.getReferenceIDIndex());

        segment.setStart((std::min)(segment.getStart(), otherSegment.getStart()));
        segment.setEnd((std::max)(segment.getEnd(), otherSegment.getEnd()));

        if (segment.getStrand() != otherSegment.getStrand()) {
            segment.setStrand(GenomicStrand::NONE);
        }
    };

    absorbSegment(sortedSegments.firstRegion, other.sortedSegments.firstRegion);
    absorbSegment(sortedSegments.secondRegion, other.sortedSegments.secondRegion);

    maxComplementarityScore = (std::max)(maxComplementarityScore, other.maxComplementarityScore);
    minHybridizationEnergy = (std::min)(minHybridizationEnergy, other.minHybridizationEnergy);

    recordIDs.insert(recordIDs.end(), other.recordIDs.begin(), other.recordIDs.end());
    complementarityScores.insert(complementarityScores.end(), other.complementarityScores.begin(),
                                 other.complementarityScores.end());
    hybridizationEnergies.insert(hybridizationEnergies.end(), other.hybridizationEnergies.begin(),
                                 other.hybridizationEnergies.end());
    crosslinkingSiteCounts.insert(crosslinkingSiteCounts.end(), other.crosslinkingSiteCounts.begin(),
                                  other.crosslinkingSiteCounts.end());

    transcriptContribution += other.transcriptContribution;
    transcriptContributionSquaredSum += other.transcriptContributionSquaredSum;
    firstArmCoverage.merge(other.firstArmCoverage);
    secondArmCoverage.merge(other.secondArmCoverage);
}

auto InteractionCluster::coverageShapeMetrics() const -> CoverageShapeMetrics {
    const ArmCoverageSummary firstSummary = firstArmCoverage.summary();
    const ArmCoverageSummary secondSummary = secondArmCoverage.summary();

    const double integratedCoverage =
        firstSummary.integratedCoverage + secondSummary.integratedCoverage;
    const double squaredCoverageIntegral =
        firstSummary.squaredCoverageIntegral + secondSummary.squaredCoverageIntegral;
    const double effectiveCoverageSpanBp =
        safeDivide(integratedCoverage * integratedCoverage, squaredCoverageIntegral);
    const size_t totalSpan = totalSpanBp();

    CoverageShapeMetrics metrics{
        .totalSpanBp = totalSpan,
        .effectiveCoverageSpanBp = effectiveCoverageSpanBp,
        .supportPerTotalBp = safeDivide(getTranscriptContribution(), static_cast<double>(totalSpan)),
        .supportPerEffectiveBp =
            safeDivide(getTranscriptContribution(), effectiveCoverageSpanBp),
        .coverageConcentration =
            clampUnitInterval(1.0 - safeDivide(effectiveCoverageSpanBp,
                                               static_cast<double>(totalSpan))),
        .coverageComponents = firstSummary.components + secondSummary.components,
        .armBalance = safeDivide(std::min(firstSummary.integratedCoverage,
                                          secondSummary.integratedCoverage),
                                 std::max(firstSummary.integratedCoverage,
                                          secondSummary.integratedCoverage)),
        .coverageProfile = {}};
    metrics.coverageProfile = classifyCoverageProfile(metrics);
    return metrics;
}

void InteractionCluster::replaceCoverageFromFragments(const RecordFragment &firstFragment,
                                                      const RecordFragment &secondFragment) {
    firstArmCoverage = {};
    secondArmCoverage = {};

    const bool firstFragmentIsFirstArm = firstFragment.genomicRegion == sortedSegments.firstRegion;
    addFragmentCoverage(firstFragmentIsFirstArm ? firstArmCoverage : secondArmCoverage,
                        firstFragment);
    addFragmentCoverage(firstFragmentIsFirstArm ? secondArmCoverage : firstArmCoverage,
                        secondFragment);
}

auto InteractionCluster::complementarityStatistics() const -> double {
    return helper::calculateMedian(complementarityScores) * maxComplementarityScore;
}

auto InteractionCluster::hybridizationEnergyStatistics() const -> double {
    return std::sqrt(helper::calculateMedian(hybridizationEnergies) * minHybridizationEnergy);
}

[[nodiscard]] auto InteractionCluster::meanCrosslinkingSiteCount() const -> double {
    return std::accumulate(crosslinkingSiteCounts.begin(), crosslinkingSiteCounts.end(), 0.0) /
           static_cast<double>(crosslinkingSiteCounts.size());
}

[[nodiscard]] auto InteractionCluster::standardDeviationCrosslinkingSiteCount() const -> double {
    if (crosslinkingSiteCounts.size() < 2) {
        return 0.0;
    }

    const double mean = meanCrosslinkingSiteCount();
    const double squaredDifferences =
        std::accumulate(crosslinkingSiteCounts.begin(), crosslinkingSiteCounts.end(), 0.0,
                        [mean](double accumulator, int32_t count) {
                            return accumulator + std::pow(count - mean, 2);
                        });
    return std::sqrt(squaredDifferences / static_cast<double>(crosslinkingSiteCounts.size()));
}

auto operator<<(std::ostream &outputStream, const InteractionCluster &interactionCluster)
    -> std::ostream & {
    outputStream << "InteractionCluster:\n"
                 << "First segment: " << interactionCluster.getFirstSegment()
                 << "Second segment: " << interactionCluster.getSecondSegment()
                 << "Count: " << std::fixed << interactionCluster.getTranscriptContribution()
                 << "\n";

    for (const auto &recordID : interactionCluster.getRecordIDs()) {
        outputStream << recordID << ", ";
    }
    outputStream << "\nComplementarity scores: ";

    for (const auto &complementarityScore : interactionCluster.getComplementarityScores()) {
        outputStream << complementarityScore << ", ";
    }
    outputStream << "\nHybridization energies: ";

    for (const auto &hybridizationEnergy : interactionCluster.getHybridizationEnergies()) {
        outputStream << hybridizationEnergy << ", ";
    }

    outputStream << "\nCrosslinking site counts: ";
    for (const auto &crosslinkingSiteCount : interactionCluster.getCrosslinkingSiteCounts()) {
        outputStream << crosslinkingSiteCount << ", ";
    }

    return outputStream;
}

}  // namespace pipelines::analyze
