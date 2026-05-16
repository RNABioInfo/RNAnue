// NOLINTBEGIN(readability-magic-numbers)

#include <gtest/gtest.h>

// Standard
#include <cmath>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Boost
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/complement.hpp>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "InteractionCluster.hpp"
#include "Region.hpp"
#include "SortedGenomicRegionPair.hpp"
#include "StatisticEvaluator.hpp"
#include "TranscriptContributionsByID.hpp"

using namespace dataTypes;
using namespace pipelines::analyze;

namespace {

auto makeCluster(std::string recordID, int start, float contribution) -> InteractionCluster {
    return InteractionCluster{
        SortedGenomicRegionPair{
            GenomicRegion{0, Region{.startPosition = start, .endPosition = start + 10},
                          GenomicStrand::FORWARD},
            GenomicRegion{1, Region{.startPosition = start + 100, .endPosition = start + 110},
                          GenomicStrand::FORWARD}},
        std::move(recordID),
        0.5,
        -10.0,
        1,
        contribution};
}

auto makeClusterWithUnitContributions(std::string recordPrefix, int start, size_t count)
    -> InteractionCluster {
    std::vector<std::string> recordIDs;
    std::vector<double> complementarityScores;
    std::vector<double> hybridizationEnergies;
    std::vector<int32_t> crosslinkingSiteCounts;

    recordIDs.reserve(count);
    complementarityScores.reserve(count);
    hybridizationEnergies.reserve(count);
    crosslinkingSiteCounts.reserve(count);

    for (size_t i = 0; i < count; ++i) {
        recordIDs.emplace_back(recordPrefix + std::to_string(i));
        complementarityScores.emplace_back(0.5);
        hybridizationEnergies.emplace_back(-10.0);
        crosslinkingSiteCounts.emplace_back(1);
    }

    return InteractionCluster{
        SortedGenomicRegionPair{
            GenomicRegion{0, Region{.startPosition = start, .endPosition = start + 10},
                          GenomicStrand::FORWARD},
            GenomicRegion{1, Region{.startPosition = start + 100, .endPosition = start + 110},
                          GenomicStrand::FORWARD}},
        std::move(recordIDs),
        std::move(complementarityScores),
        std::move(hybridizationEnergies),
        std::move(crosslinkingSiteCounts),
        static_cast<float>(count)};
}

auto makeAnnotated(InteractionCluster cluster, std::string firstFeatureID,
                   std::string secondFeatureID) -> AnnotatedInteractionCluster {
    return AnnotatedInteractionCluster{
        std::move(cluster), std::move(firstFeatureID), std::move(secondFeatureID)};
}

auto binomialUpperTail(double trials, double probability, double observedEffectiveSuccesses)
    -> double {
    const auto distribution = boost::math::binomial_distribution<double>(trials, probability);
    return boost::math::cdf(boost::math::complement(distribution,
                                                    observedEffectiveSuccesses - 1.0));
}

}  // namespace

TEST(StatisticEvaluatorTests, IntegerContributionsUseAbundanceCorrectedBinomialTail) {
    std::vector<AnnotatedInteractionCluster> clusters{
        makeAnnotated(makeClusterWithUnitContributions("ab", 0, 8), "A", "B"),
        makeAnnotated(makeClusterWithUnitContributions("cc", 20, 2), "C", "C"),
    };
    const TranscriptContributionsByID background{{"A", 9.5F}, {"B", 9.5F}, {"C", 9.5F}};

    const auto results = StatisticEvaluator::evaluate(clusters, background, 1.0);

    ASSERT_EQ(results.size(), 2UL);
    EXPECT_NEAR(results[0].getPValue(), binomialUpperTail(10.0, 2.0 / 3.0, 8.0), 1e-12);
    EXPECT_NEAR(results[1].getPValue(), binomialUpperTail(10.0, 1.0 / 3.0, 2.0), 1e-12);
}

TEST(StatisticEvaluatorTests, FractionalContributionsUseEffectiveTrialCount) {
    InteractionCluster fractionalCluster = makeCluster("frac1", 0, 0.5F);
    EXPECT_TRUE(fractionalCluster.merge(makeCluster("frac2", 0, 0.5F),
                                        GenomicStrandSpecificity::SPECIFIC));

    std::vector<AnnotatedInteractionCluster> clusters{
        makeAnnotated(std::move(fractionalCluster), "A", "B"),
        makeAnnotated(makeCluster("unit", 20, 1.0F), "C", "C"),
    };
    const TranscriptContributionsByID background{{"A", 9.5F}, {"B", 9.5F}, {"C", 9.5F}};

    const auto results = StatisticEvaluator::evaluate(clusters, background, 1.0);

    const double totalContribution = 2.0;
    const double squaredContributionSum = 1.5;
    const double effectiveTrials = totalContribution * totalContribution / squaredContributionSum;
    const double effectiveSuccessScale = totalContribution / squaredContributionSum;

    ASSERT_EQ(results.size(), 2UL);
    EXPECT_NEAR(results[0].getPValue(),
                binomialUpperTail(effectiveTrials, 2.0 / 3.0, 1.0 * effectiveSuccessScale),
                1e-12);
    EXPECT_NEAR(results[1].getPValue(),
                binomialUpperTail(effectiveTrials, 1.0 / 3.0, 1.0 * effectiveSuccessScale),
                1e-12);
}

TEST(StatisticEvaluatorTests, MultipleClustersForSameFeaturePairSplitPairWeight) {
    std::vector<AnnotatedInteractionCluster> clusters{
        makeAnnotated(makeCluster("ab1", 0, 1.0F), "A", "B"),
        makeAnnotated(makeCluster("ab2", 20, 1.0F), "A", "B"),
        makeAnnotated(makeCluster("cc", 40, 1.0F), "C", "C"),
    };
    const TranscriptContributionsByID background{{"A", 9.5F}, {"B", 9.5F}, {"C", 9.5F}};

    const auto results = StatisticEvaluator::evaluate(clusters, background, 1.0);

    ASSERT_EQ(results.size(), 3UL);
    for (const auto &result : results) {
        EXPECT_NEAR(result.getPValue(), binomialUpperTail(3.0, 1.0 / 3.0, 1.0), 1e-12);
    }
}

TEST(StatisticEvaluatorTests, BackgroundSmoothingKeepsMissingFeaturesTestable) {
    std::vector<AnnotatedInteractionCluster> clusters{
        makeAnnotated(makeCluster("missing1", 0, 1.0F), "missingA", "missingB"),
        makeAnnotated(makeCluster("missing2", 20, 1.0F), "missingC", "missingC"),
    };
    const TranscriptContributionsByID background{};

    const auto results = StatisticEvaluator::evaluate(clusters, background, 1.0);

    ASSERT_EQ(results.size(), 2UL);
    for (const auto &result : results) {
        EXPECT_TRUE(std::isfinite(result.getPValue()));
        EXPECT_TRUE(std::isfinite(result.getPadj()));
    }
}

TEST(StatisticEvaluatorTests, BenjaminiHochbergThresholdFiltersAdjustedPValues) {
    std::vector<AnnotatedInteractionCluster> clusters{
        makeAnnotated(makeClusterWithUnitContributions("ab", 0, 8), "A", "B"),
        makeAnnotated(makeClusterWithUnitContributions("cc", 20, 2), "C", "C"),
    };
    const TranscriptContributionsByID background{{"A", 9.5F}, {"B", 9.5F}, {"C", 9.5F}};

    const auto results = StatisticEvaluator::evaluate(clusters, background, 0.7);

    ASSERT_EQ(results.size(), 1UL);
    EXPECT_EQ(results.front().getFirstFeatureID(), "A");
    EXPECT_EQ(results.front().getSecondFeatureID(), "B");
    EXPECT_NEAR(results.front().getPadj(), 2.0 * binomialUpperTail(10.0, 2.0 / 3.0, 8.0),
                1e-12);
}

TEST(StatisticEvaluatorTests, TranscriptContributionParserReadsAllRows) {
    std::istringstream input{"A\t1.5\nB\t2\nA\t0.5\n"};
    TranscriptContributionsByID transcriptCounts;

    parseTranscriptContributionStream(input, transcriptCounts);

    EXPECT_FLOAT_EQ(transcriptCounts.at("A"), 2.0F);
    EXPECT_FLOAT_EQ(transcriptCounts.at("B"), 2.0F);
}

// NOLINTEND(readability-magic-numbers)
