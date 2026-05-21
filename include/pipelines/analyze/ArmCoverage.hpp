#pragma once

// Standard
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace pipelines::analyze {

struct CoverageInterval {
    int32_t start{};
    int32_t end{};

    auto operator==(const CoverageInterval&) const -> bool = default;
};

struct CoverageRun {
    int32_t start{};
    int32_t end{};
    double coverage{};

    auto operator==(const CoverageRun&) const -> bool = default;
};

struct ArmCoverageSummary {
    double integratedCoverage{};
    double squaredCoverageIntegral{};
    double maxCoverage{};
    size_t components{};

    [[nodiscard]] auto effectiveSpan() const noexcept -> double {
        if (integratedCoverage <= 0.0 || squaredCoverageIntegral <= 0.0) {
            return 0.0;
        }

        return integratedCoverage * integratedCoverage / squaredCoverageIntegral;
    }
};

class ArmCoverage {
   public:
    void addInterval(int32_t start, int32_t end, double weight) {
        if (start >= end || !std::isfinite(weight) || weight <= 0.0) {
            return;
        }

        events.push_back({.position = start, .delta = weight});
        events.push_back({.position = end, .delta = -weight});
    }

    void addIntervals(const std::vector<CoverageInterval>& intervals, double weight) {
        for (const auto& interval : intervals) {
            addInterval(interval.start, interval.end, weight);
        }
    }

    void merge(const ArmCoverage& other) {
        events.insert(events.end(), other.events.begin(), other.events.end());
    }

    [[nodiscard]] auto runs() const -> std::vector<CoverageRun> {
        if (events.empty()) {
            return {};
        }

        std::vector<Event> sortedEvents = events;
        std::ranges::sort(sortedEvents, [](const Event& lhs, const Event& rhs) {
            return lhs.position < rhs.position;
        });

        std::vector<CoverageRun> coverageRuns;
        coverageRuns.reserve(sortedEvents.size());

        double currentCoverage = 0.0;
        int32_t previousPosition = sortedEvents.front().position;
        size_t index = 0;

        while (index < sortedEvents.size()) {
            const int32_t position = sortedEvents[index].position;

            if (position > previousPosition && currentCoverage > coverageEpsilon) {
                addRun(coverageRuns, previousPosition, position, currentCoverage);
            }

            double delta = 0.0;
            while (index < sortedEvents.size() && sortedEvents[index].position == position) {
                delta += sortedEvents[index].delta;
                ++index;
            }

            currentCoverage += delta;
            if (std::abs(currentCoverage) <= coverageEpsilon) {
                currentCoverage = 0.0;
            }

            previousPosition = position;
        }

        return coverageRuns;
    }

    [[nodiscard]] auto summary() const -> ArmCoverageSummary {
        const auto coverageRuns = runs();
        return summarizeRuns(coverageRuns);
    }

    [[nodiscard]] static auto summarizeRuns(const std::vector<CoverageRun>& coverageRuns)
        -> ArmCoverageSummary {
        ArmCoverageSummary summary{};

        for (const auto& run : coverageRuns) {
            if (run.start >= run.end || run.coverage <= coverageEpsilon) {
                continue;
            }

            const auto length = static_cast<double>(run.end - run.start);
            summary.integratedCoverage += length * run.coverage;
            summary.squaredCoverageIntegral += length * run.coverage * run.coverage;
            summary.maxCoverage = std::max(summary.maxCoverage, run.coverage);
        }

        const double componentThreshold = componentCoverageThreshold(summary.maxCoverage);
        bool inComponent = false;
        int32_t previousEnd = 0;

        for (const auto& run : coverageRuns) {
            if (run.start >= run.end) {
                continue;
            }

            const bool aboveThreshold = run.coverage >= componentThreshold;
            if (!aboveThreshold) {
                inComponent = false;
                previousEnd = run.end;
                continue;
            }

            if (!inComponent || run.start > previousEnd) {
                ++summary.components;
            }

            inComponent = true;
            previousEnd = run.end;
        }

        return summary;
    }

   private:
    struct Event {
        int32_t position{};
        double delta{};
    };

    static constexpr double coverageEpsilon = 1e-12;

    std::vector<Event> events;

    static void addRun(std::vector<CoverageRun>& coverageRuns, int32_t start, int32_t end,
                       double coverage) {
        if (!coverageRuns.empty() && coverageRuns.back().end == start &&
            std::abs(coverageRuns.back().coverage - coverage) <= coverageEpsilon) {
            coverageRuns.back().end = end;
            return;
        }

        coverageRuns.push_back({.start = start, .end = end, .coverage = coverage});
    }

    [[nodiscard]] static auto componentCoverageThreshold(double maxCoverage) noexcept -> double {
        return std::max(1.0, 0.05 * maxCoverage);
    }
};

}  // namespace pipelines::analyze
