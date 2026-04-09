#pragma once

// Standard
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <ranges>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

// seqan3
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

// Internal
#include "FigurePlotter.hpp"
#include "HitGroupFailureReason.hpp"
#include "ReadGroup.hpp"
#include "ReadGroupPreprocessor.hpp"
#include "SamRecord.hpp"
#include "SegemehlHitGroupConstructor.hpp"

namespace pipelines::detect {

namespace fs = std::filesystem;

struct SegemehlReadGroupPreprocessorMetrics {
    [[nodiscard]] constexpr auto getPassedHitGroupsPerReadGroup() const noexcept
        -> const std::vector<size_t>& {
        return passedHitGroupsPerReadGroup;
    }
    [[nodiscard]] constexpr auto getFailedHitGroupsPerReadGroup() const noexcept
        -> const std::vector<size_t>& {
        return failedHitGroupsPerReadGroup;
    }
    [[nodiscard]] constexpr auto getSuccesReadGroupCount() const -> size_t {
        return passedReadGroupCount;
    }
    [[nodiscard]] constexpr auto getTotalFailedReadGroupCount() const -> size_t {
        return totalFailReadGroupCount;
    }
    [[nodiscard]] constexpr auto getTotalReadGroupCount() const -> size_t {
        return passedReadGroupCount + totalFailReadGroupCount;
    }

    void operator+=(const SegemehlReadGroupPreprocessorMetrics& other) noexcept {
        passedHitGroupsPerReadGroup.insert(passedHitGroupsPerReadGroup.end(),
                                           other.passedHitGroupsPerReadGroup.begin(),
                                           other.passedHitGroupsPerReadGroup.end());
        failedHitGroupsPerReadGroup.insert(failedHitGroupsPerReadGroup.end(),
                                           other.failedHitGroupsPerReadGroup.begin(),
                                           other.failedHitGroupsPerReadGroup.end());

        passedReadGroupCount += other.passedReadGroupCount;
        totalFailReadGroupCount += other.totalFailReadGroupCount;
    }

    void createPlots(const fs::path& outDir) const {
        auto plotter = plotting::FigurePlotter(
            plotting::FigureConfig::makeDefault("Pre Processing Hit Group Contribution", {outDir}));

        plotter.addStackedBar<size_t>(
            {.title = "Failed / Passed Read Groups", .xlabel = "Sample", .ylabel = "Count"},
            {.data = {{passedReadGroupCount}, {totalFailReadGroupCount}},
             .legendTitle = "Status",
             .legendLabels = {"Passed", "Failed"},
             .groupLabels = std::nullopt});

        plotter.addScatter<size_t>({.title = "Passed vs. Failed Hit Groups per Read Group",
                                    .xlabel = "Count Failed Hit Group",
                                    .ylabel = "Count Passed Hit Group"},
                                   {.x_vals = failedHitGroupsPerReadGroup,
                                    .y_vals = passedHitGroupsPerReadGroup,
                                    .datapointLabel = "Read Group"});

        plotter.save();
    }

    void addContext(ConstructedEvaluationContextVariant& contextVariant) noexcept {
        std::visit(
            [&](const auto& context) {
                using CtxT = std::remove_cvref_t<decltype(context)>;
                if constexpr (CtxT::isFailed()) {
                    ++currentFailedHitGroupCount;
                } else {
                    ++currentPassedHitGroupCount;
                }
            },
            contextVariant);
    }

    void finalizeReadGroup() noexcept {
        passedHitGroupsPerReadGroup.push_back(currentPassedHitGroupCount);
        failedHitGroupsPerReadGroup.push_back(currentFailedHitGroupCount);

        // Checks if we had at least one passed hit group in read group
        if (currentPassedHitGroupCount == 0) {
            ++totalFailReadGroupCount;
        } else {
            ++passedReadGroupCount;
        }

        currentPassedHitGroupCount = 0;
        currentFailedHitGroupCount = 0;
    }

    void finalizeReadGroupAllFailed() noexcept {
        size_t sumHitGroups = currentPassedHitGroupCount + currentFailedHitGroupCount;

        passedHitGroupsPerReadGroup.push_back(0);
        failedHitGroupsPerReadGroup.push_back(sumHitGroups);

        ++totalFailReadGroupCount;

        currentPassedHitGroupCount = 0;
        currentFailedHitGroupCount = 0;
    }

   private:
    std::vector<size_t> passedHitGroupsPerReadGroup;
    std::vector<size_t> failedHitGroupsPerReadGroup;
    size_t passedReadGroupCount{0};
    size_t totalFailReadGroupCount{0};

    size_t currentPassedHitGroupCount{0};
    size_t currentFailedHitGroupCount{0};
};

struct SegemehlReadGroupPreprocessorConfig {
    SegemehlHitGroupConstructionParameters constructionParams;
    size_t maxPrimaryAlignmentCount;
};

struct SegemehlReadGroupPreprocessor {
    SegemehlReadGroupPreprocessor(SegemehlReadGroupPreprocessorConfig config)
        : collationConfig(config) {}

    [[nodiscard]] constexpr auto getMetrics() const -> const SegemehlReadGroupPreprocessorMetrics& {
        return metrics;
    }

    auto operator()(ReadGroup&& readGroup) -> std::vector<ConstructedEvaluationContextVariant> {
        std::unordered_map<int, SegemehlHitGroupConstructor> hitGroupConstructorByTag;
        hitGroupConstructorByTag.reserve(readGroup.size());

        for (auto&& record : std::move(readGroup)) {
            assert(record.tags().contains("HI"_tag));
            const int hitGroupTag = record.tags().get<"HI"_tag>();

            auto [iterator, inserted] = hitGroupConstructorByTag.try_emplace(
                hitGroupTag,
                SegemehlHitGroupConstructor{std::move(record), collationConfig.constructionParams});

            // If it already existed, just insert the moved record
            if (!inserted) {
                iterator->second.insert(std::move(record));
            }
        }

        std::vector<ConstructedEvaluationContextVariant> constructedContexts;
        constructedContexts.reserve(hitGroupConstructorByTag.size());
        for (auto& [tag, constructor] : hitGroupConstructorByTag) {
            constructedContexts.push_back(constructor.getConstructedEvalContext());

            metrics.addContext(constructedContexts.back());
        }

        if (primaryAlignmentCount(constructedContexts) <=
            collationConfig.maxPrimaryAlignmentCount) {
            metrics.finalizeReadGroup();
            return constructedContexts;
        }

        // In case there are more alignments then the cutoff allows for
        std::vector<ConstructedEvaluationContextVariant> errorContextVariants;
        errorContextVariants.reserve(constructedContexts.size());

        for (auto&& contextVariant : constructedContexts) {
            std::visit(
                [&](auto&& context) {
                    using ctx_t = decltype(context);

                    if constexpr (!context.isFailed()) {
                        errorContextVariants.push_back(
                            std::forward<ctx_t>(context).with(HitGroupFailureReason::MULTIMAPPING));
                    } else {
                        errorContextVariants.push_back(std::forward<ctx_t>(context));
                    }
                },
                std::move(contextVariant));
        }

        return errorContextVariants;
    }

   private:
    SegemehlReadGroupPreprocessorConfig collationConfig;
    SegemehlReadGroupPreprocessorMetrics metrics{};

    [[nodiscard]] static auto primaryAlignmentCount(
        const std::vector<ConstructedEvaluationContextVariant>& contextVariants) -> size_t {
        std::unordered_map<size_t, size_t> freq;
        freq.reserve(contextVariants.size());
        size_t maxCount = 0;

        for (auto const& contextVariant : contextVariants) {
            size_t mapQ = std::visit(
                [](auto const& ctx) {
                    return ctx.group->getRecordContainer().front().mapping_quality();
                },
                contextVariant);

            size_t& cnt = freq[mapQ];
            ++cnt;
            maxCount = std::max(cnt, maxCount);
        }

        assert(maxCount > 0);
        return maxCount;
    }
};

static_assert(ReadGroupPreprocessor<SegemehlReadGroupPreprocessor>);

}  // namespace pipelines::detect
