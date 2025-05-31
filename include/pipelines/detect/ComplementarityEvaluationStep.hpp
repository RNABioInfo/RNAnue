#pragma once

// Standard
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

// seqan3
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

// matplot
#include <matplot/matplot.h>

// Internal
#include "CoOptimalPairwiseAligner.hpp"
#include "EvaluationContext.hpp"
#include "HitGroup.hpp"
#include "ReadGroupEvaluationParameters.hpp"
#include "SamRecord.hpp"
#include "SplitRecords.hpp"

namespace pipelines::detect {

using namespace dataTypes;
namespace fs = std::filesystem;

struct ComplementarityEvaluationResult {
    bool passed;
    std::optional<CoOptimalPairwiseAligner::Result> alignmentResult;

    void addTags(SamRecord &record) const {
        if (!alignmentResult) {
            return;
        }

        const auto &result = *alignmentResult;

        const int length = static_cast<int>(result.endPositions.first) -
                           static_cast<int>(result.beginPositions.first);
        record.tags()["XL"_tag] = length;
        record.tags()["XC"_tag] = static_cast<float>(result.complementarity);
        record.tags()["XR"_tag] = static_cast<float>(result.fraction);
        record.tags()["XS"_tag] = result.score;
    }
};

struct ComplementarityEvaluationStepConfig {
    double minComplementarity;
    double minFraction;

    template <ReadGroupEvaluationParameters::Type EvalParamT>
    [[nodiscard]] static auto makeConfig(const EvalParamT &params)
        -> ComplementarityEvaluationStepConfig {
        if constexpr (std::same_as<EvalParamT, ReadGroupEvaluationParameters::Base>) {
            return {.minComplementarity = params.minComplementarity,
                    .minFraction = params.minComplementarityFraction};
        } else {
            return {.minComplementarity = params.baseParameters.minComplementarity,
                    .minFraction = params.baseParameters.minComplementarityFraction};
        }
    }
};

struct ComplementarityEvaluationStepMetrics {
   public:
    explicit ComplementarityEvaluationStepMetrics(ComplementarityEvaluationStepConfig config)
        : config(config) {}

    [[nodiscard]] constexpr auto getPassedCount() const noexcept -> size_t { return passedCount; }
    [[nodiscard]] constexpr auto getFailedCount() const noexcept -> size_t { return failedCount; }
    [[nodiscard]] constexpr auto getComplementarityScores() const noexcept
        -> const std::vector<float> & {
        return complementarityScores;
    }
    [[nodiscard]] constexpr auto getFractionScores() const noexcept -> const std::vector<float> & {
        return fractionScores;
    }

    void operator+=(const ComplementarityEvaluationResult &result) noexcept;
    void operator+=(const ComplementarityEvaluationStepMetrics &other) noexcept;
    void createPlots(const fs::path &outDir) const;

   private:
    ComplementarityEvaluationStepConfig config;

    size_t passedCount{0}, failedCount{0};
    std::vector<float> complementarityScores;
    std::vector<float> fractionScores;
};

class ComplementarityEvaluationStep {
   public:
    ComplementarityEvaluationStep(ComplementarityEvaluationStepConfig config)
        : config(config), metrics(config) {}

    using step_result_t = ComplementarityEvaluationResult;

    template <typename... OtherResults>
    using step_return_t = std::variant<
        std::vector<EvaluationContext<SingletonHitGroup, OtherResults..., step_result_t>>,
        EvaluationContext<ChimericHitGroup, OtherResults..., step_result_t>>;

    template <typename... OtherResults>
    auto operator()(EvaluationContext<ChimericHitGroup, OtherResults...> &&ctx)
        -> step_return_t<OtherResults...> {
        auto result = evaluate(ctx.group->getRecordContainer());
        metrics += result;

        if (!result.passed) {
            return std::forward_like<decltype(ctx)>(ctx).with(result).getDeconstructedContexts();
        }

        return std::forward_like<decltype(ctx)>(ctx).with(result);
    };

    [[nodiscard]] auto getMetrics() const -> const ComplementarityEvaluationStepMetrics & {
        return metrics;
    }

   private:
    ComplementarityEvaluationStepConfig config;
    ComplementarityEvaluationStepMetrics metrics;

    [[nodiscard]] auto evaluate(const ChimericRecords &splitRecords) const noexcept
        -> ComplementarityEvaluationResult;

    [[nodiscard]] auto isPassingFilters(
        const CoOptimalPairwiseAligner::Result &alignmentResult) const noexcept -> bool;

    static constexpr auto complementaryScoringScheme() -> seqan3::nucleotide_scoring_scheme<int8_t>;
};

}  // namespace pipelines::detect
