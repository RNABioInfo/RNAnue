#pragma once

// Standard
#include <concepts>
#include <cstddef>
#include <memory>
#include <tuple>
#include <type_traits>
#include <variant>
#include <vector>

// Internal
#include "HitGroup.hpp"
#include "HitGroupFailureReason.hpp"
#include "SamRecord.hpp"
#include "TupleForEach.hpp"

namespace pipelines::detect {

using namespace dataTypes;

template <typename Group, typename... EvaluationResultsT>
struct EvaluationContext {
    std::unique_ptr<Group> group;
    std::tuple<EvaluationResultsT...> evaluationResults;

    using hit_group_t = Group;

    template <typename T>
    static constexpr bool has_tag_v = (std::is_same_v<T, EvaluationResultsT> || ...);

    explicit constexpr EvaluationContext(std::unique_ptr<Group> group) noexcept
        requires(sizeof...(EvaluationResultsT) == 0)
        : group(std::move(group)) {}

    constexpr EvaluationContext(std::unique_ptr<Group> group,
                                std::tuple<EvaluationResultsT...>&& evalResults) noexcept
        : group(std::move(group)), evaluationResults(std::move(evalResults)) {}

    constexpr EvaluationContext(std::unique_ptr<Group> group,
                                const std::tuple<EvaluationResultsT...>& evalResults) noexcept
        : group(std::move(group)), evaluationResults(evalResults) {}

    EvaluationContext(EvaluationContext const& other)
        // if other.group is non-null, copy-construct *other.group, otherwise leave nullptr
        : group(other.group ? std::make_unique<Group>(*other.group) : nullptr),
          evaluationResults(other.evaluationResults)  // tuple copy
    {}

    [[nodiscard]] auto getRecords() const noexcept -> Group::RecordsContainerT& {
        auto& records = group->getRecordContainer();
        const auto& results = evaluationResults;
        size_t index = 0;

        for (SamRecord& record : records) {
            tupleForEach(results, [&](const auto& res) {
                if constexpr (requires { res.addTags(record); }) {
                    res.addTags(record);
                } else if constexpr (requires { res.addTags(record, index); }) {
                    res.addTags(record, index);
                } else {
                    static_assert("Step result does not implement addTags");
                }
            });

            ++index;
        }
        return records;
    }

    template <class U>
        requires std::movable<std::decay_t<U>>
    auto with(U&& newValue) && -> EvaluationContext<Group, EvaluationResultsT..., std::decay_t<U>> {
        using T = std::decay_t<U>;
        auto newExtras =
            std::tuple_cat(std::move(evaluationResults), std::tuple<T>(std::forward<U>(newValue)));
        return EvaluationContext<Group, EvaluationResultsT..., T>{std::move(group),
                                                                  std::move(newExtras)};
    }

    template <typename T>
    [[nodiscard]] constexpr auto hasTag() const noexcept -> bool {
        return has_tag_v<T>;
    }

    template <std::size_t I>
    constexpr auto get() & noexcept -> decltype(auto) {
        return std::get<I>(evaluationResults);
    }

    template <std::size_t I>
    constexpr auto get() const& noexcept -> decltype(auto) {
        return std::get<I>(evaluationResults);
    }

    template <typename T>
    constexpr auto get() & noexcept -> decltype(auto) {
        return std::get<T>(evaluationResults);
    }

    template <typename T>
    constexpr auto get() const& noexcept -> decltype(auto) {
        return std::get<T>(evaluationResults);
    }

    static consteval auto isFailed() -> bool { return has_tag_v<HitGroupFailureReason>; }

    [[nodiscard]] auto getDeconstructedContexts() const noexcept
        -> std::vector<EvaluationContext<SingletonHitGroup, EvaluationResultsT...>>
        requires(!std::same_as<hit_group_t, SingletonHitGroup>)
    {
        std::vector<EvaluationContext<SingletonHitGroup, EvaluationResultsT...>> contexts;
        contexts.reserve(group->getRecordCount());

        for (auto&& newGroup : group->getDeconstructedHitGroups()) {
            contexts.emplace_back(
                std::make_unique<SingletonHitGroup>(std::forward<decltype(newGroup)>(newGroup)),
                evaluationResults);
        }

        return contexts;
    }
};

using EvaluationContextVariant =
    std::variant<EvaluationContext<SingletonHitGroup>, EvaluationContext<ChimericHitGroup>,
                 EvaluationContext<MultimericHitGroup>>;

using ErrorEvaluationContextVariant =
    std::variant<EvaluationContext<SingletonHitGroup, HitGroupFailureReason>,
                 EvaluationContext<ChimericHitGroup, HitGroupFailureReason>,
                 EvaluationContext<MultimericHitGroup, HitGroupFailureReason>>;

}  // namespace pipelines::detect
