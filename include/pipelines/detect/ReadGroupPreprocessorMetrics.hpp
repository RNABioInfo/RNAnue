#pragma once
#include <map>

// Standard
#include <cassert>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <type_traits>
#include <variant>
#include <vector>

// seqan3
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

// Internal
#include "FigurePlotter.hpp"
#include "HitGroupFailureReason.hpp"
#include "ReadGroupPreprocessor.hpp"

namespace pipelines::detect {
namespace fs = std::filesystem;
struct ReadGroupPreprocessorMetrics {
    std::map<dataTypes::HitGroupFailureReason, size_t> failures;

    [[nodiscard]] constexpr auto getPassedHitGroupsPerReadGroup() const
        -> const std::vector<size_t>& {
        return passedHitGroupsPerReadGroup;
    }
    [[nodiscard]] constexpr auto getFailedHitGroupsPerReadGroup() const
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

    void operator+=(const ReadGroupPreprocessorMetrics& other) {
        for (const auto& [reason, count] : other.failures) failures[reason] += count;
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

    void addContext(ConstructedEvaluationContextVariant& contextVariant) {
        std::visit(
            [&](const auto& context) {
                using CtxT = std::remove_cvref_t<decltype(context)>;
                if constexpr (CtxT::isFailed()) {
                    ++currentFailedHitGroupCount;
                    ++failures[context.template get<dataTypes::HitGroupFailureReason>()];
                } else {
                    ++currentPassedHitGroupCount;
                }
            },
            contextVariant);
    }

    void finalizeReadGroup() {
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

    void finalizeReadGroupAllFailed() {
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

}  // namespace pipelines::detect
