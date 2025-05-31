#include "ComplementarityEvaluationStep.hpp"

// Standard
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <optional>
#include <ranges>
#include <tuple>
#include <vector>

// seqan3
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

// Internal
#include "CoOptimalPairwiseAligner.hpp"
#include "FigurePlotter.hpp"
#include "FileFormat.hpp"
#include "SplitRecords.hpp"
#include "UnderlyingSequence.hpp"

namespace pipelines::detect {

void ComplementarityEvaluationStepMetrics::operator+=(
    const ComplementarityEvaluationResult &result) noexcept {
    if (result.passed) {
        passedCount++;
    } else {
        failedCount++;
    }

    if (!result.alignmentResult) {
        return;
    }

    complementarityScores.push_back(result.alignmentResult->complementarity);
    fractionScores.push_back(result.alignmentResult->fraction);
}

void ComplementarityEvaluationStepMetrics::operator+=(
    const ComplementarityEvaluationStepMetrics &other) noexcept {
    passedCount += other.passedCount;
    failedCount += other.failedCount;

    complementarityScores.insert(complementarityScores.end(), other.complementarityScores.begin(),
                                 other.complementarityScores.end());
    fractionScores.insert(fractionScores.end(), other.fractionScores.begin(),
                          other.fractionScores.end());
}

void ComplementarityEvaluationStepMetrics::createPlots(const fs::path &outDir) const {
    auto plotter = plotting::FigurePlotter(plotting::FigureConfig::makeDefault(
        "Complementarity Metrics of all Chimeric Hits", {outDir}));

    // Failed Passed Plot
    std::vector<std::vector<size_t>> failedPassed = {{failedCount}, {passedCount}};
    plotter.addStackedBar<size_t>(
        {.title = "Failed / Passed Hit Groups", .xlabel = "Sample", .ylabel = "Count"},
        {.data = failedPassed,
         .legendTitle = "Evalutation Result",
         .legendLabels = {"Passed", "Failed"},
         .groupLabels = std::nullopt});

    // Complementarity vs. Fraction Plot
    const std::vector<double> x_vals = {complementarityScores.begin(), complementarityScores.end()};
    const std::vector<double> y_vals = {fractionScores.begin(), fractionScores.end()};
    plotter.addScatter<double>({.title = "Alignment Complementarity vs. Fraction Matched",
                                .xlabel = "Fraction of Matched vs. Alignment Length",
                                .ylabel = "Fraction of Matched vs. Shortest Record Length"},
                               {.x_vals = x_vals, .y_vals = y_vals, .datapointLabel = "Hit Group"});

    // Complementarity Plot
    plotter.addHistogram<float>({.title = "Alignment Complementarity Scores",
                                 .xlabel = "Fraction of Matched vs. Alignment Length",
                                 .ylabel = "Count"},
                                plotting::HistogramData{.data = complementarityScores,
                                                        .datapointLabel = "Hit Groups",
                                                        .cutoffLabel = "Min: ",
                                                        .cutoffValue = config.minComplementarity});

    // Fraction Plot
    plotter.addHistogram<float>({.title = "Fraction Matched Scores",
                                 .xlabel = "Fraction of Matched vs. Shortest Record Length",
                                 .ylabel = "Count"},
                                plotting::HistogramData{.data = fractionScores,
                                                        .datapointLabel = "Hit Groups",
                                                        .cutoffLabel = "Min: ",
                                                        .cutoffValue = config.minFraction});

    plotter.save();
}

auto ComplementarityEvaluationStep::evaluate(const ChimericRecords &splitRecords) const noexcept
    -> ComplementarityEvaluationResult {
    const auto sequence1View = splitRecords.first().sequence() | std::views::reverse |
                               views::underlying_sequence(splitRecords.first().flag());
    const auto sequence2View =
        splitRecords.second().sequence() | views::underlying_sequence(splitRecords.second().flag());

    const auto sequence1 = std::vector<seqan3::dna5>(sequence1View.begin(), sequence1View.end());
    const auto sequence2 = std::vector<seqan3::dna5>(sequence2View.begin(), sequence2View.end());

    const auto results =
        CoOptimalPairwiseAligner{ComplementarityEvaluationStep::complementaryScoringScheme()}
            .getLocalAlignments(std::tie(sequence1, sequence2));

    if (results.empty()) {
        return {.passed = false, .alignmentResult = std::nullopt};
    }

    auto alignmentResult = *std::ranges::max_element(results, std::less{});

    return {.passed = isPassingFilters(alignmentResult), .alignmentResult = alignmentResult};
}

constexpr auto ComplementarityEvaluationStep::complementaryScoringScheme()
    -> seqan3::nucleotide_scoring_scheme<int8_t> {
    using namespace seqan3::literals;

    seqan3::nucleotide_scoring_scheme scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}};

    scheme.score('A'_dna5, 'T'_dna5) = 1;
    scheme.score('T'_dna5, 'A'_dna5) = 1;
    scheme.score('G'_dna5, 'C'_dna5) = 1;
    scheme.score('C'_dna5, 'G'_dna5) = 1;
    scheme.score('G'_dna5, 'T'_dna5) = 1;
    scheme.score('T'_dna5, 'G'_dna5) = 1;

    scheme.score('T'_dna5, 'T'_dna5) = -1;
    scheme.score('A'_dna5, 'A'_dna5) = -1;
    scheme.score('C'_dna5, 'C'_dna5) = -1;
    scheme.score('G'_dna5, 'G'_dna5) = -1;

    scheme.score('A'_dna5, 'G'_dna5) = -1;
    scheme.score('A'_dna5, 'C'_dna5) = -1;
    scheme.score('G'_dna5, 'A'_dna5) = -1;
    scheme.score('C'_dna5, 'A'_dna5) = -1;

    scheme.score('T'_dna5, 'C'_dna5) = -1;
    scheme.score('C'_dna5, 'T'_dna5) = -1;

    return scheme;
}

auto ComplementarityEvaluationStep::isPassingFilters(
    const CoOptimalPairwiseAligner::Result &alignmentResult) const noexcept -> bool {
    return alignmentResult.complementarity >= config.minComplementarity &&
           alignmentResult.fraction >= config.minFraction;
};

}  // namespace pipelines::detect
