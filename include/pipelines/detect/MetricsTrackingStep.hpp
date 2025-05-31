#pragma once

#include <concepts>
#include <cstddef>
#include <ostream>
#include <string>
#include <unordered_map>

#include "EvaluationContext.hpp"
#include "HitGroup.hpp"
#include "SplitRecords.hpp"
namespace pipelines::detect {

using namespace dataTypes;

struct Metrics {
    using CountsByRecordType = std::unordered_map<SplitRecordType, size_t>;

    [[nodiscard]] auto getReadCount() const noexcept -> size_t { return readCount; }
    [[nodiscard]] auto getHitGroupCounts() const noexcept -> const CountsByRecordType& {
        return hitGroupCount;
    }

    void incrementHitGroups(const SplitRecordType& splitRecordType) {
        hitGroupCount[splitRecordType] += 1;
    }
    void incrementReads() noexcept { readCount++; }

    void operator+=(const Metrics& other) {
        readCount += other.getReadCount();
        for (const auto& [recordType, count] : other.getHitGroupCounts()) {
            hitGroupCount[recordType] += count;
        }
    }

   private:
    CountsByRecordType hitGroupCount;
    size_t readCount{0};
};

struct MetricsTrackingStepResult {
    void addTags(auto& /* unused */) const noexcept {}
};

struct MetricsTrackingStep {
    using step_result_t = MetricsTrackingStepResult;

    template <typename HitGroupT, typename... OtherResultsT>
    auto operator()(EvaluationContext<HitGroupT, OtherResultsT...>&& ctx) {
        using ctx_t = decltype(ctx);

        metrics.incrementHitGroups(HitGroupT::hitGroupSplitRecordType);

        if (currentID != ctx.group->readID()) {
            currentID = ctx.group->readID();
            metrics.incrementReads();
        }

        return std::forward_like<ctx_t>(ctx).with(MetricsTrackingStepResult{});
    }

    [[nodiscard]] auto getMetrics() const noexcept -> const Metrics& { return metrics; }

    void operator+=(const MetricsTrackingStep& other) { metrics += other.getMetrics(); }

   private:
    std::string currentID;

    Metrics metrics{};
};

inline auto operator<<(std::ostream& out, const MetricsTrackingStep& step) -> std::ostream& {
    out << "Read count (" << step.getMetrics().getReadCount() << "); ";

    size_t totalCount = 0;

    for (const auto& [type, count] : step.getMetrics().getHitGroupCounts()) {
        out << type << " hit group count (" << count << "); ";

        totalCount += count;
    }

    out << "Total hit group count (" << totalCount << ")";

    return out;
}

}  // namespace pipelines::detect
