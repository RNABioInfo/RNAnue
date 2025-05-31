#pragma once

// Standard
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <optional>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <variant>
#include <vector>

// ViennaRNA
extern "C" {
#include <ViennaRNA/cofold.h>
#include <ViennaRNA/subopt.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
}

// Internal
#include "EvaluationContext.hpp"
#include "FigurePlotter.hpp"
#include "HitGroup.hpp"
#include "HitGroupEvaluationResult.hpp"
#include "ReadGroupEvaluationParameters.hpp"
#include "SamRecord.hpp"
#include "SplitRecords.hpp"

namespace pipelines::detect {
namespace fs = std::filesystem;
using namespace seqan3::literals;
using namespace dataTypes;

struct HybridizationEvaluationStepConfig {
    double mfeThreshold;
    bool includeWobbleBasePairsInCrosslinkingSites;

    template <ReadGroupEvaluationParameters::Type EvalParamT>
    [[nodiscard]] static auto makeConfig(const EvalParamT &params)
        -> HybridizationEvaluationStepConfig {
        if constexpr (std::same_as<EvalParamT, ReadGroupEvaluationParameters::Base>) {
            return HybridizationEvaluationStepConfig{
                params.mfeThreshold, params.includeWobbleBasePairsInCrosslinkingSites};
        } else {
            return HybridizationEvaluationStepConfig{
                params.baseParameters.mfeThreshold,
                params.baseParameters.includeWobbleBasePairsInCrosslinkingSites};
        }
    }
};

struct HybridizationEvaluationResult {
    bool passed;
    std::optional<double> energy;
    std::optional<HitGroupEvaluation::CrosslinkingResult> crosslinkingResult;

    void addTags(SamRecord &record, size_t index) const {
        if (energy) {
            record.tags()["XE"_tag] = static_cast<float>(*energy);
        }

        if (crosslinkingResult) {
            record.tags()["XD"_tag] = crosslinkingResult->getDotbracketString();
            record.tags()["XO"_tag] =
                static_cast<int32_t>(crosslinkingResult->getInterCrosslinkingCount());

            record.tags()["XA"_tag] = crosslinkingResult->getIntraSequenceCrosslinking(index);

            record.tags()["XI"_tag] = crosslinkingResult->getInterSequenceCrosslinking(index);
        }
    }
};

struct HybridizationEvaluationStepMetrics {
   public:
    HybridizationEvaluationStepMetrics(HybridizationEvaluationStepConfig config) : config(config) {}

    [[nodiscard]] constexpr auto getPassedCount() const noexcept -> size_t { return passedCount; }
    [[nodiscard]] constexpr auto getFailedCount() const noexcept -> size_t { return failedCount; }
    [[nodiscard]] constexpr auto getMFEValues() const noexcept -> const std::vector<double> & {
        return mfeValues;
    }

    template <typename T>
    [[nodiscard]] constexpr auto getIntraCrosslinkingCounts() const noexcept -> std::vector<T> {
        std::vector<T> counts;
        // convert each char to size_t
        std::ranges::transform(crosslinkingResults, std::back_inserter(counts),
                               [](const HitGroupEvaluation::CrosslinkingResult &result) -> T {
                                   return static_cast<T>(result.getIntraCrosslinkingCount());
                               });

        return counts;
    }

    template <typename T>
    [[nodiscard]] constexpr auto getInterCrosslinkingCounts() const noexcept -> std::vector<T> {
        std::vector<T> counts;
        // convert each char to size_t
        std::ranges::transform(crosslinkingResults, std::back_inserter(counts),
                               [](const HitGroupEvaluation::CrosslinkingResult &result) -> T {
                                   return static_cast<T>(result.getInterCrosslinkingCount());
                               });

        return counts;
    }

    void operator+=(const HybridizationEvaluationResult &result) noexcept {
        if (result.passed) {
            passedCount++;
        } else {
            failedCount++;
        }

        if (result.energy) {
            mfeValues.push_back(*result.energy);
        }

        if (result.crosslinkingResult) {
            crosslinkingResults.push_back(*result.crosslinkingResult);
        }
    };

    void operator+=(const HybridizationEvaluationStepMetrics &other) noexcept {
        passedCount += other.passedCount;
        failedCount += other.failedCount;

        mfeValues.insert(mfeValues.end(), other.mfeValues.begin(), other.mfeValues.end());
        crosslinkingResults.insert(crosslinkingResults.end(), other.crosslinkingResults.begin(),
                                   other.crosslinkingResults.end());
    };

    void createPlots(const fs::path &outDir) const {
        auto plotter = plotting::FigurePlotter(plotting::FigureConfig::makeDefault(
            "Hybridization Metrics of all Chimeric Hits", outDir));

        // Failed Passed Plot
        std::vector<std::vector<size_t>> failedPassed = {{failedCount}, {passedCount}};
        plotter.addStackedBar<size_t>(
            {.title = "Failed / Passed Hit Groups", .xlabel = "Sample", .ylabel = "Count"},
            {.data = failedPassed,
             .legendTitle = "Evalutation Result",
             .legendLabels = {"Passed", "Failed"},
             .groupLabels = std::nullopt});

        const std::vector<double> y_vals = getInterCrosslinkingCounts<double>();
        plotter.addScatter<double>(
            {.title = "Intermolecular Crosslinks vs. MFE",
             .xlabel = "Minimum Free Energy",
             .ylabel = "Intermolecular Crosslink Count"},
            {.x_vals = mfeValues, .y_vals = y_vals, .datapointLabel = "Hit Group"});

        plotter.addHistogram<double>({.title = "Minimum Free Energy Scores",
                                      .xlabel = "Minimum Free Energy",
                                      .ylabel = "Count"},
                                     plotting::HistogramData{.data = mfeValues,
                                                             .datapointLabel = "Hit Groups",
                                                             .cutoffLabel = "Max: ",
                                                             .cutoffValue = config.mfeThreshold});

        plotter.addHistogram<size_t>(
            {.title = "Predicted Crosslinking Sites",
             .xlabel = "Crosslinking Count",
             .ylabel = "Count"},
            plotting::HistogramData{.data = getInterCrosslinkingCounts<size_t>(),
                                    .datapointLabel = "Intermolecular",
                                    .cutoffLabel = std::nullopt,
                                    .cutoffValue = std::nullopt},
            plotting::HistogramData{.data = getIntraCrosslinkingCounts<size_t>(),
                                    .datapointLabel = "Intramolecular",
                                    .cutoffLabel = std::nullopt,
                                    .cutoffValue = std::nullopt});

        plotter.save();
    };

   private:
    HybridizationEvaluationStepConfig config;

    size_t passedCount{0}, failedCount{0};
    std::vector<double> mfeValues;
    std::vector<HitGroupEvaluation::CrosslinkingResult> crosslinkingResults;
};

class HybridizationEvaluationStep {
   public:
    HybridizationEvaluationStep(HybridizationEvaluationStepConfig config)
        : config(config), metrics(config) {};

    using step_result_t = HybridizationEvaluationResult;
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

    [[nodiscard]] auto getMetrics() const -> const HybridizationEvaluationStepMetrics & {
        return metrics;
    }

   private:
    [[nodiscard]] auto evaluate(const ChimericRecords &splitRecords) const
        -> HybridizationEvaluationResult;
    HybridizationEvaluationStepConfig config;
    HybridizationEvaluationStepMetrics metrics;

    [[nodiscard]] auto isPassingFilters(double energy) const noexcept -> bool;
};

}  // namespace pipelines::detect
