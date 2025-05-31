#pragma once

// Internal
#include <cmath>
#include <cstdlib>
#include <optional>
#include <utility>
#include <variant>

#include "HitGroup.hpp"
#include "HitGroupEvaluationResult.hpp"
#include "HitGroupFilterReason.hpp"
#include "SplitRecords.hpp"

namespace dataTypes::HitGroupEvaluation {

template <RecordsContainer T>
class EvaluatedHitGroup : public HitGroup<T> {
   public:
    // The constructor requires the underlying HitGroup and both evaluation results.
    EvaluatedHitGroup(HitGroup<T>&& hitGroup, ComplementarityResult&& compResult,
                      HybridizationResult&& hybResult)
        : HitGroup<T>(std::move(hitGroup)),
          complementarityResult(std::move(compResult)),
          hybridizationResult(std::move(hybResult)) {}

    EvaluatedHitGroup(HitGroup<T>&& hitGroup, ComplementarityResult&& compResult,
                      HybridizationResult&& hybResult, const HitGroupFilterReason filterReason)
        : HitGroup<T>(std::move(hitGroup)),
          complementarityResult(std::move(compResult)),
          hybridizationResult(std::move(hybResult)) {
        this->filterReason = filterReason;
    }

    // For instance, an operator< that uses evaluation results might be provided.
    // (You can adjust the details as needed; here we illustrate using a threshold.)
    auto operator<(const EvaluatedHitGroup& other) const -> bool {
        constexpr double DIFF_THRESHOLD = 0.05;

        if (std::abs(complementarityResult.fraction - other.complementarityResult.fraction) >=
            DIFF_THRESHOLD) {
            return complementarityResult.fraction < other.complementarityResult.fraction;
        }
        if (std::abs(complementarityResult.complementarity -
                     other.complementarityResult.complementarity) >= DIFF_THRESHOLD) {
            return complementarityResult.complementarity <
                   other.complementarityResult.complementarity;
        }
        if (std::abs(hybridizationResult.energy - other.hybridizationResult.energy) >=
            DIFF_THRESHOLD) {
            return hybridizationResult.energy < other.hybridizationResult.energy;
        }
        if (hybridizationResult.crosslinkingResult.has_value() &&
            other.hybridizationResult.crosslinkingResult.has_value()) {
            return hybridizationResult.crosslinkingResult->getTotalCrosslinkingCount() <
                   other.hybridizationResult.crosslinkingResult->getTotalCrosslinkingCount();
        }
        return false;
    }

    auto operator>(const EvaluatedHitGroup& other) const -> bool { return other < *this; }

    [[nodiscard]] constexpr auto isValid() const noexcept -> bool { !filterReason.has_value(); }

    [[nodiscard]] constexpr auto getFilterReason() const noexcept
        -> std::optional<HitGroupFilterReason> {
        return filterReason;
    }

   private:
    ComplementarityResult complementarityResult;
    HybridizationResult hybridizationResult;

    std::optional<HitGroupFilterReason> filterReason;
};

// First, for the case when no extra evaluation is required.
// This is the case for SingletonRecord.
template <>
class EvaluatedHitGroup<SingletonRecord> : public HitGroup<SingletonRecord> {
   public:
    using Base = HitGroup<SingletonRecord>;

    // Inherit the constructors from the base hitgroup.
    using Base::HitGroup;

    [[nodiscard]] static constexpr auto isValid() noexcept -> bool { return true; }
};

template <>
class EvaluatedHitGroup<MultimericRecords> : public HitGroup<MultimericRecords> {
   public:
    using Base = HitGroup<MultimericRecords>;

    // Inherit the constructors from the base hitgroup.
    using Base::HitGroup;

    // TODO: Currently multisplits are not supported
    [[nodiscard]] static constexpr auto isValid() noexcept -> bool { return false; }
};

using EvaluatedHitGroupVariant =
    std::variant<EvaluatedHitGroup<SingletonRecord>, EvaluatedHitGroup<ChimericRecords>,
                 EvaluatedHitGroup<MultimericRecords>>;

}  // namespace dataTypes::HitGroupEvaluation
