// NOLINTBEGIN(readability-magic-numbers)

#include <gtest/gtest.h>

// Standard
#include <string>
#include <vector>

// Internal
#include "ArmCoverage.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "InteractionCluster.hpp"
#include "ParseSamRecords.hpp"
#include "RecordFragment.hpp"
#include "Region.hpp"

using namespace dataTypes;
using namespace pipelines::analyze;

namespace {

auto makeFragment(std::string recordID, int referenceID, int start, int end,
                  float contribution = 1.0F) -> RecordFragment {
    return {.recordID = std::move(recordID),
            .genomicRegion =
                GenomicRegion{referenceID, Region{.startPosition = start, .endPosition = end},
                              GenomicStrand::FORWARD},
            .complementarityScore = 0.5,
            .hybridizationEnergy = -1.0,
            .interCrosslinkingSiteCount = 1,
            .transcriptContribution = contribution,
            .coverageIntervals = {{.start = start, .end = end}}};
}

auto makeCluster(std::string recordID, int firstStart, int firstEnd, int secondStart,
                 int secondEnd, float contribution = 1.0F) -> InteractionCluster {
    return InteractionCluster::fromRecordFragments(
        makeFragment(recordID, 0, firstStart, firstEnd, contribution),
        makeFragment(recordID, 1, secondStart, secondEnd, contribution));
}

}  // namespace

TEST(ArmCoverageTests, WeightedRunsAndEffectiveSpanReflectCoverageShape) {
    ArmCoverage coverage;
    coverage.addInterval(0, 10, 1.0);
    coverage.addInterval(5, 15, 1.0);

    const auto runs = coverage.runs();
    ASSERT_EQ(runs.size(), 3UL);
    EXPECT_EQ(runs[0].start, 0);
    EXPECT_EQ(runs[0].end, 5);
    EXPECT_DOUBLE_EQ(runs[0].coverage, 1.0);
    EXPECT_EQ(runs[1].start, 5);
    EXPECT_EQ(runs[1].end, 10);
    EXPECT_DOUBLE_EQ(runs[1].coverage, 2.0);
    EXPECT_EQ(runs[2].start, 10);
    EXPECT_EQ(runs[2].end, 15);
    EXPECT_DOUBLE_EQ(runs[2].coverage, 1.0);

    const auto summary = coverage.summary();
    EXPECT_DOUBLE_EQ(summary.integratedCoverage, 20.0);
    EXPECT_DOUBLE_EQ(summary.squaredCoverageIntegral, 30.0);
    EXPECT_NEAR(summary.effectiveSpan(), 400.0 / 30.0, 1e-12);
    EXPECT_DOUBLE_EQ(summary.maxCoverage, 2.0);
    EXPECT_EQ(summary.components, 1UL);
}

TEST(ArmCoverageTests, CigarAwareFragmentCoverageSkipsDeletionsAndReferenceSkips) {
    const auto rawSam = R"(@HD	VN:1.6
@SQ	SN:ref	LN:1000
read1	0	ref	1	255	5M2D5M3N5M	*	0	0	ACGTACGTACGTACG	*	XE:f:-1	XC:f:0.5	XO:i:1	XB:f:0.5
)";

    const auto records = parseSamRecords(rawSam);
    ASSERT_EQ(records.size(), 1UL);

    const auto fragment = RecordFragment::fromSamRecord(records.front());
    ASSERT_TRUE(fragment.has_value());

    ASSERT_EQ(fragment->coverageIntervals.size(), 3UL);
    EXPECT_EQ(fragment->coverageIntervals[0], (CoverageInterval{.start = 0, .end = 5}));
    EXPECT_EQ(fragment->coverageIntervals[1], (CoverageInterval{.start = 7, .end = 12}));
    EXPECT_EQ(fragment->coverageIntervals[2], (CoverageInterval{.start = 15, .end = 20}));
}

TEST(ArmCoverageTests, EffectiveSupportDensityIsStableAcrossGappedMerges) {
    auto mergedCluster = makeCluster("left", 0, 10, 1000, 1010);
    const auto distantCluster = makeCluster("right", 100, 110, 1100, 1110);

    mergedCluster.absorbValidatedComponentMember(distantCluster);

    const auto metrics = mergedCluster.coverageShapeMetrics();
    EXPECT_EQ(metrics.totalSpanBp, 220UL);
    EXPECT_NEAR(metrics.effectiveCoverageSpanBp, 40.0, 1e-12);
    EXPECT_NEAR(metrics.supportPerTotalBp, 2.0 / 220.0, 1e-12);
    EXPECT_NEAR(metrics.supportPerEffectiveBp, 2.0 / 40.0, 1e-12);
    EXPECT_EQ(metrics.coverageComponents, 4UL);
    EXPECT_EQ(metrics.coverageProfile, "multi_peak_refine");
}

TEST(ArmCoverageTests, CoverageProfilesClassifyExpectedShapes) {
    EXPECT_EQ(makeCluster("compact", 0, 10, 100, 110)
                  .coverageShapeMetrics()
                  .coverageProfile,
              "compact_dense");

    EXPECT_EQ(makeCluster("broad", 0, 200, 1000, 1200, 20.0F)
                  .coverageShapeMetrics()
                  .coverageProfile,
              "broad_diffuse");

    EXPECT_EQ(makeCluster("low", 0, 200, 1000, 1200)
                  .coverageShapeMetrics()
                  .coverageProfile,
              "low_support_density");

    EXPECT_EQ(makeCluster("imbalanced", 0, 100, 1000, 1010, 10.0F)
                  .coverageShapeMetrics()
                  .coverageProfile,
              "imbalanced_support");

    auto multiPeak = makeCluster("left", 0, 10, 1000, 1010);
    multiPeak.absorbValidatedComponentMember(makeCluster("right", 100, 110, 1100, 1110));
    EXPECT_EQ(multiPeak.coverageShapeMetrics().coverageProfile, "multi_peak_refine");
}

// NOLINTEND(readability-magic-numbers)
