#pragma once

// Standard
#include <optional>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>

// Internal
#include "CooptimalPairwiseAligner.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"

namespace pipelines::detect {

class SplitRecordsComplementarityEvaluator {
   public:
    using Result = CoOptimalPairwiseAligner::Result;

    SplitRecordsComplementarityEvaluator() = delete;
    ~SplitRecordsComplementarityEvaluator() = delete;
    SplitRecordsComplementarityEvaluator(const SplitRecordsComplementarityEvaluator &) = delete;
    auto operator=(const SplitRecordsComplementarityEvaluator &)
        -> SplitRecordsComplementarityEvaluator & = delete;
    SplitRecordsComplementarityEvaluator(SplitRecordsComplementarityEvaluator &&) = delete;
    auto operator=(SplitRecordsComplementarityEvaluator &&)
        -> SplitRecordsComplementarityEvaluator & = delete;

    static auto evaluate(const SplitRecords &splitRecords,
                         const SplitRecordsEvaluationParameters::BaseParameters &parameters)
        -> std::optional<SplitRecordsComplementarityEvaluator::Result>;

   private:
    static constexpr auto complementaryScoringScheme() -> seqan3::nucleotide_scoring_scheme<int8_t>;

    static auto getOptimalAlignment(
        const std::vector<CoOptimalPairwiseAligner::Result> &alignResults,
        double minComplementarity, double minComplementarityFraction)
        -> std::optional<SplitRecordsComplementarityEvaluator::Result>;
};

}  // namespace pipelines::detect
