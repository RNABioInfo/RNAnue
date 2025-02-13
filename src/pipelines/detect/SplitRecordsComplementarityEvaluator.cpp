#include "SplitRecordsComplementarityEvaluator.hpp"

// Standard
#include <cstdint>
#include <optional>
#include <ranges>
#include <tuple>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

// Internal
#include "CoOptimalPairwiseAligner.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"
#include "UnderlyingSequence.hpp"

namespace pipelines::detect {

auto SplitRecordsComplementarityEvaluator::evaluate(
    const SplitRecords &splitRecords,
    const SplitRecordsEvaluationParameters::BaseParameters &parameters)
    -> std::optional<SplitRecordsComplementarityEvaluator::Result> {
    const auto &record1 = splitRecords[0];
    const auto &record2 = splitRecords[1];

    const auto sequence1View =
        record1.sequence() | std::views::reverse | views::underlying_sequence(record1.flag());
    const auto sequence2View = record2.sequence() | views::underlying_sequence(record1.flag());

    const auto sequence1 = std::vector<seqan3::dna5>(sequence1View.begin(), sequence1View.end());
    const auto sequence2 = std::vector<seqan3::dna5>(sequence2View.begin(), sequence2View.end());

    const auto results =
        CoOptimalPairwiseAligner{SplitRecordsComplementarityEvaluator::complementaryScoringScheme()}
            .getLocalAlignments(std::tie(sequence1, sequence2));

    return SplitRecordsComplementarityEvaluator::getOptimalAlignment(
        results, parameters.minComplementarity, parameters.minComplementarityFraction);
}

constexpr auto SplitRecordsComplementarityEvaluator::complementaryScoringScheme()
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

auto SplitRecordsComplementarityEvaluator::getOptimalAlignment(
    const std::vector<CoOptimalPairwiseAligner::Result> &alignResults,
    const double minComplementarity, const double minComplementarityFraction)
    -> std::optional<SplitRecordsComplementarityEvaluator::Result> {
    if (alignResults.empty()) {
        return std::nullopt;
    }

    std::optional<CoOptimalPairwiseAligner::Result> optimalAlignment = std::nullopt;

    for (const auto &alignResult : alignResults) {
        if (alignResult.complementarity < minComplementarity ||
            alignResult.fraction < minComplementarityFraction) {
            continue;
        }

        if (!optimalAlignment.has_value() ||
            alignResult.complementarity > optimalAlignment.value().complementarity ||
            (alignResult.complementarity == optimalAlignment.value().complementarity &&
             alignResult.fraction > optimalAlignment.value().fraction)) {
            optimalAlignment = alignResult;
        }
    }

    return optimalAlignment;
}

}  // namespace pipelines::detect
