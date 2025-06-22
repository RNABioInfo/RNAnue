#pragma once

// Standard
#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <iterator>
#include <numeric>
#include <optional>
#include <ranges>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

// seqan3
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

// Internal
#include "EvaluationContext.hpp"
#include "FigurePlotter.hpp"
#include "SamRecord.hpp"
#include "SplitRecords.hpp"
#include "Utility.hpp"
#include "VariantMeta.hpp"

namespace pipelines::detect {

namespace fs = std::filesystem;
using namespace dataTypes;

template <typename C>
concept HasRecordContainer = requires(C const& ctx) {
    { ctx.group->getRecordContainer() } -> std::ranges::input_range;
};

struct ReadGroupPostScoringResult {
    double contributionScore;

    void addTags(SamRecord& record) const noexcept {
        record.tags()["XB"_tag] = static_cast<float>(contributionScore);
    }
};

struct ReadGroupPostScoringStepConfig {
    double minContribution;
};

struct ReadGroupPostScoringStepMetrics {
   public:
    explicit ReadGroupPostScoringStepMetrics(ReadGroupPostScoringStepConfig config)
        : config(config) {}

    [[nodiscard]] constexpr auto getHitGroupsPreFilter() const noexcept
        -> const std::vector<size_t>& {
        return hitGroupsPreFilter;
    }
    [[nodiscard]] constexpr auto getHitGroupsPostFilter() const noexcept
        -> const std::vector<size_t>& {
        return hitGroupsPostFilter;
    }
    [[nodiscard]] constexpr auto getHitGroupContributionScores() const noexcept
        -> const std::vector<double>& {
        return hitGroupContributionScores;
    }

    [[nodiscard]] constexpr auto getContributionScoreSum() const noexcept -> double {
        return std::reduce(hitGroupContributionScores.begin(), hitGroupContributionScores.end());
    }

    [[nodiscard]] constexpr auto getReadGroupCount() const noexcept -> size_t {
        return readGroupCount;
    }

    [[nodiscard]] constexpr auto getContributionScoreByRecordType() const noexcept
        -> const std::unordered_map<SplitRecordType, double>& {
        return contributionScoreByRecordType;
    }

    void incrementReadGroupCount() noexcept { ++readGroupCount; }

    void addPreFilterCount(const size_t count) noexcept { hitGroupsPreFilter.push_back(count); }

    void addPostFilterCount(const size_t count) noexcept { hitGroupsPostFilter.push_back(count); }

    void addHitGroupContributionScore(const double score) noexcept {
        hitGroupContributionScores.push_back(score);
    }

    template <typename... EvalContextT>
    void addHitGroupCounts(const std::vector<std::variant<EvalContextT...>>& contexts) {
        for (const auto& contextVariant : contexts) {
            std::visit(
                [&](const auto& context) {
                    using CtxT = std::remove_cvref_t<decltype(context)>;

                    static constexpr SplitRecordType splitRecordType =
                        CtxT::hit_group_t::hitGroupSplitRecordType;

                    const ReadGroupPostScoringResult& result =
                        context.template get<ReadGroupPostScoringResult>();

                    if (!std::isnan(result.contributionScore)) {
                        contributionScoreByRecordType[splitRecordType] += result.contributionScore;
                    }
                },
                contextVariant);
        }
    }

    void operator+=(const ReadGroupPostScoringStepMetrics& other) noexcept {
        hitGroupsPreFilter.insert(hitGroupsPreFilter.end(), other.hitGroupsPreFilter.begin(),
                                  other.hitGroupsPreFilter.end());
        hitGroupsPostFilter.insert(hitGroupsPostFilter.end(), other.hitGroupsPostFilter.begin(),
                                   other.hitGroupsPostFilter.end());
        hitGroupContributionScores.insert(hitGroupContributionScores.end(),
                                          other.hitGroupContributionScores.begin(),
                                          other.hitGroupContributionScores.end());

        readGroupCount += other.readGroupCount;

        for (const auto& [typeIndex, count] : other.contributionScoreByRecordType) {
            contributionScoreByRecordType[typeIndex] += count;
        }
    }

    void createPlots(const fs::path& outDir) const {
        auto plotter = plotting::FigurePlotter(
            plotting::FigureConfig::makeDefault("Post Processing Hit Group Contributions", outDir));

        std::vector<std::string> labels;
        std::vector<std::vector<double>> contributionScores;
        std::ranges::for_each(contributionScoreByRecordType, [&](const auto& map) {
            labels.insert(labels.begin(), to_string(map.first));
            contributionScores.push_back({static_cast<double>(map.second)});
        });

        plotter.addStackedBar<double>({.title = "Contributions by Record Type",
                                       .xlabel = "Sample",
                                       .ylabel = "Contributions"},
                                      {.data = {contributionScores},
                                       .legendTitle = "Evalutation Result",
                                       .legendLabels = labels,
                                       .groupLabels = std::nullopt});

        plotter.addHistogram<size_t>({.title = "Hit Groups / Read Pre Contribution Filter",
                                      .xlabel = "Hit Groups In Read",
                                      .ylabel = "Count"},
                                     plotting::HistogramData{.data = hitGroupsPreFilter,
                                                             .datapointLabel = "Hit Groups",
                                                             .cutoffLabel = std::nullopt,
                                                             .cutoffValue = std::nullopt});

        plotter.addHistogram<size_t>({.title = "Hit Groups / Read Post Contribution Filter",
                                      .xlabel = "Hit Groups In Read",
                                      .ylabel = "Count"},
                                     plotting::HistogramData{.data = hitGroupsPostFilter,
                                                             .datapointLabel = "Hit Groups",
                                                             .cutoffLabel = std::nullopt,
                                                             .cutoffValue = std::nullopt});

        plotter.addHistogram<double>(
            {.title = "Hit Groups Contribution Scores",
             .xlabel = "Contribution Score",
             .ylabel = "Count"},
            plotting::HistogramData{.data = hitGroupContributionScores,
                                    .datapointLabel = "Hit Group",
                                    .cutoffLabel = "Min: ",
                                    .cutoffValue = config.minContribution});

        plotter.save();
    }

   private:
    ReadGroupPostScoringStepConfig config;

    std::vector<size_t> hitGroupsPreFilter;
    std::vector<size_t> hitGroupsPostFilter;
    std::vector<double> hitGroupContributionScores;
    std::unordered_map<SplitRecordType, double> contributionScoreByRecordType;
    size_t readGroupCount{0};
};

struct ReadGroupPostScoringStep {
    explicit ReadGroupPostScoringStep(ReadGroupPostScoringStepConfig config)
        : config(config), metrics(config) {}

    template <typename... ContextVariants>
        requires(... && HasRecordContainer<ContextVariants>)
    [[nodiscard]] auto operator()(std::vector<std::variant<ContextVariants...>>&& contexts) {
        using result_t = std::variant<std::decay_t<decltype(std::declval<ContextVariants>().with(
            ReadGroupPostScoringResult{}))>...>;

        if (contexts.empty()) {
            return std::vector<result_t>{};
        }

        metrics.addPreFilterCount(contexts.size());

        const std::vector<double> editScores = getEditDistanceScores(getEditDistances(contexts));
        const std::vector<double> alignScores = getAlignScores(getAlignedLengths(contexts));

        const std::vector<IndexedScore> finalScores = getFinalScores(editScores, alignScores);

        metrics.addPostFilterCount(finalScores.size());

        // Only if after scoring there are contexts available the read group is deemed passing
        if (!finalScores.empty()) {
            metrics.incrementReadGroupCount();
        }

        std::vector<result_t> newContexts;
        newContexts.reserve(finalScores.size());

        for (const auto& [index, score] : finalScores) {
            std::visit(
                [&](auto&& ctx) {
                    newContexts.emplace_back(
                        std::forward<decltype(ctx)>(ctx).with(ReadGroupPostScoringResult{score}));
                },
                std::move(contexts[index]));
        }

        // TODO: Assign directly inside new context assignment
        metrics.addHitGroupCounts(newContexts);

        return newContexts;
    };

    [[nodiscard]] auto getMetrics() const -> const ReadGroupPostScoringStepMetrics& {
        return metrics;
    }

   private:
    ReadGroupPostScoringStepConfig config;
    ReadGroupPostScoringStepMetrics metrics;

    struct IndexedScore {
        size_t index;
        double score;
    };

    template <typename Group, typename... OtherResults>
    auto transform(EvaluationContext<Group, OtherResults...>&& ctx)
        -> EvaluationContext<Group, OtherResults..., ReadGroupPostScoringResult>;

    template <typename... ContextVariants>
        requires(... && HasRecordContainer<ContextVariants>)
    [[nodiscard]] static auto getEditDistances(
        const std::vector<std::variant<ContextVariants...>>& contexts) noexcept
        -> std::vector<size_t> {
        std::vector<size_t> editDistances;
        editDistances.reserve(contexts.size());

        std::ranges::transform(
            contexts, std::back_inserter(editDistances),
            [&](const std::variant<ContextVariants...>& contextVariant) -> size_t {
                return std::visit(
                    [&](const auto& context) -> size_t { return getEditDistance(context); },
                    contextVariant);
            },
            std::identity{});

        return editDistances;
    }

    [[nodiscard]] static auto getEditDistance(const auto& context) noexcept -> size_t {
        const auto& container = context.group->getRecordContainer();

        return std::transform_reduce(
            container.begin(), container.end(), std::size_t{0}, std::plus<>{},
            [](SamRecord const& record) { return record.tags().template get<"NM"_tag>(); });
    }

    template <typename... ContextVariants>
        requires(... && HasRecordContainer<ContextVariants>)
    [[nodiscard]] static auto getAlignedLengths(
        const std::vector<std::variant<ContextVariants...>>& contexts) noexcept
        -> std::vector<size_t> {
        std::vector<size_t> alignmentLengths;
        alignmentLengths.reserve(contexts.size());

        std::ranges::transform(
            contexts, std::back_inserter(alignmentLengths),
            [&](const std::variant<ContextVariants...>& contextVariant) -> size_t {
                return std::visit(
                    [&](const auto& context) -> size_t { return getAlignedLength(context); },
                    contextVariant);
            },
            std::identity{});

        return alignmentLengths;
    }

    [[nodiscard]] static auto getAlignedLength(const auto& context) noexcept -> size_t {
        const auto& container = context.group->getRecordContainer();

        return std::transform_reduce(
            container.begin(), container.end(), std::size_t{0}, std::plus<>{},
            [](SamRecord const& record) { return alignmentLength(record); });
    }

    [[nodiscard]] static auto getEditDistanceScores(const std::vector<size_t>& editDistances)
        -> std::vector<double> {
        if (editDistances.size() == 1) {
            return {1};
        }

        const size_t sumEditDistances = std::reduce(editDistances.begin(), editDistances.end());
        const size_t denomEditDistances = sumEditDistances + editDistances.size();

        std::vector<double> editDistanceScores;
        editDistanceScores.reserve(editDistances.size());

        std::ranges::transform(editDistances, std::back_inserter(editDistanceScores),
                               [&](const size_t editDistance) -> double {
                                   return 1 - ((static_cast<double>(editDistance) + 1) /
                                               static_cast<double>(denomEditDistances));
                               });

        return editDistanceScores;
    }

    [[nodiscard]] static auto getAlignScores(const std::vector<size_t>& alignmentLengths)
        -> std::vector<double> {
        if (alignmentLengths.size() == 1) {
            return {1};
        }

        const size_t sumAlignLengths =
            std::reduce(alignmentLengths.begin(), alignmentLengths.end());

        std::vector<double> alignScores;
        alignScores.reserve(alignmentLengths.size());

        std::ranges::transform(
            alignmentLengths, std::back_inserter(alignScores),
            [&](const size_t alignLength) -> double {
                return (static_cast<double>(alignLength) / static_cast<double>(sumAlignLengths));
            });

        return alignScores;
    }

    [[nodiscard]] auto getFinalScores(const std::same_as<std::vector<double>> auto&... scores)
        -> std::vector<IndexedScore> {
        std::vector<IndexedScore> indexedScores;
        double sumOfValidScores{0};

        size_t index{0};

        for (const auto& elems : std::ranges::zip_view(scores...)) {
            // First, check if any score in the tuple is NaN. If so, skip this group.
            bool anyNan =
                std::apply([](auto... vals) { return ((std::isnan(vals)) || ...); }, elems);
            if (anyNan) {
                ++index;
                continue;
            }

            double totalScore = std::apply([](auto... vals) { return (vals + ...); }, elems);
            double rescaledScore = totalScore / sizeof...(scores);

            // Ensure the score meets the minimum and isn't NaN (should never be due to our earlier
            // check)
            if (rescaledScore < config.minContribution || std::isnan(rescaledScore)) {
                ++index;
                continue;
            }

            sumOfValidScores += rescaledScore;
            indexedScores.emplace_back(index, rescaledScore);

            ++index;
        }

        if (helper::isApproxEqual(0.0, sumOfValidScores) || sumOfValidScores < 0.0) {
            return {};
        }

        // Rescale all scores to sum up to 1
        double rescalingFactor = 1 / sumOfValidScores;

        for (auto& indexedScore : indexedScores) {
            indexedScore.score *= rescalingFactor;

            metrics.addHitGroupContributionScore(indexedScore.score);
        }

        return indexedScores;
    }
};

// static_assert(ReadGroupPostprocessor<ReadGroupPostScoringStep, EvaluationContextVariant>);

}  // namespace pipelines::detect
