#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include <matplot/axes_objects/histogram.h>

#include "HistogramBins.hpp"

namespace {
void expectValidEdges(const std::vector<double>& values) {
    const auto edges = plotting::histogramBinEdges(values);
    ASSERT_GE(edges.size(), 2);
    EXPECT_LE(edges.size(), 65537);
    EXPECT_LE(edges.front(), *std::ranges::min_element(values));
    EXPECT_GE(edges.back(), *std::ranges::max_element(values));
    for (size_t index = 1; index < edges.size(); ++index) {
        EXPECT_TRUE(std::isfinite(edges[index - 1]));
        EXPECT_GT(edges[index], edges[index - 1]);
    }
    EXPECT_TRUE(std::isfinite(edges.back()));
    const auto counts = matplot::histogram::histogram_count(values, edges);
    EXPECT_EQ(std::accumulate(counts.begin(), counts.end(), size_t{0}), values.size());
}
}  // namespace

TEST(HistogramBins, EmptyObservationsHaveNoEdges) {
    EXPECT_TRUE(plotting::histogramBinEdges({}).empty());
}

TEST(HistogramBins, SingletonAndConstantObservationsUseOneNonzeroWidthBin) {
    for (const auto& values : std::vector<std::vector<double>>{
             {0}, {1}, {-8.75}, {0, 0, 0}, {1, 1, 1, 1}}) {
        SCOPED_TRACE(::testing::PrintToString(values));
        expectValidEdges(values);
        EXPECT_EQ(plotting::histogramBinEdges(values).size(), 2);
    }
}

TEST(HistogramBins, SmallAndSkewedSamplesKeepAllObservationsInValidBins) {
    for (const auto& values : std::vector<std::vector<double>>{
             {0, 1}, {0, 1, 2}, {0, 0, 0, 1}, {-10.5, -10, -8.75, -6.25},
             {1, std::nextafter(1.0, 2.0)},
             {-std::numeric_limits<double>::max(), std::numeric_limits<double>::max()},
             {std::numeric_limits<double>::max()}}) {
        SCOPED_TRACE(::testing::PrintToString(values));
        expectValidEdges(values);
    }
}

TEST(HistogramBins, NonFiniteObservationsAreRejectedBeforeBinning) {
    EXPECT_THROW(plotting::histogramBinEdges({std::numeric_limits<double>::quiet_NaN()}),
                 std::invalid_argument);
    EXPECT_THROW(plotting::histogramBinEdges({std::numeric_limits<double>::infinity()}),
                 std::invalid_argument);
}
