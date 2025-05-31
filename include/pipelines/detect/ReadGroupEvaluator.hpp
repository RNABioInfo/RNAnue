#pragma once

// Standard
#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <iterator>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

// seqan3
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

// Internal
#include "EvaluationContext.hpp"
#include "HitGroup.hpp"
#include "HitGroupFailureReason.hpp"
#include "IsVariant.hpp"
#include "ReadGroup.hpp"
#include "ReadGroupPreprocessor.hpp"
#include "VariantMeta.hpp"
#include "VariantUnion.hpp"

namespace pipelines::detect {

struct IdentityStep {
    using step_result_t = HitGroupFailureReason;

    template <typename G, typename... Es>
    auto operator()(EvaluationContext<G, Es...>&& ctx) const noexcept {
        auto new_ctx =
            std::forward<decltype(ctx)>(ctx).with(HitGroupFailureReason::FAILED_COMPLEMENTARITY);
        return new_ctx;
    }
};

template <typename EvalStepT, typename HitGroupT, typename... EvalResultT>
concept EvalStep =
    requires(EvalStepT evalStep, EvaluationContext<HitGroupT, EvalResultT...>&& evalContext) {
        {
            std::move(evalStep)(std::move(evalContext))
        } noexcept -> std::same_as<
            EvaluationContext<HitGroupT, EvalResultT..., typename EvalStepT::step_result_t>>;
        typename EvalStepT::step_result_t;
    };

static_assert(EvalStep<IdentityStep, SingletonHitGroup>);
static_assert(EvalStep<IdentityStep, ChimericHitGroup>);
static_assert(EvalStep<IdentityStep, MultimericHitGroup>);

template <typename EvalContextVariantT, typename FailedEvalContextT>
struct EvalResult {
    std::vector<EvalContextVariantT> successes;
    std::vector<FailedEvalContextT> failures;
};

template <ReadGroupPreprocessor Preprocessor, typename Postprocessor, typename... EvalSteps>
class ReadGroupEvaluator {
   private:
    // Deduced types for Evaluation bases on EvalSteps
    using HitGroupTypeList =
        detail::meta::list<SingletonHitGroup, ChimericHitGroup, MultimericHitGroup>;

    using FailedCollationResultT = detail::meta::contexts_variant_t<
        HitGroupTypeList, detail::meta::list<detail::meta::list<HitGroupFailureReason>>>;

    using SuccessCollationResultT =
        detail::meta::contexts_variant_t<HitGroupTypeList,
                                         detail::meta::list<detail::meta::list<>>>;

    using EvaluationStepResultVariantT =
        detail::meta::apply_all_steps_t<SuccessCollationResultT, EvalSteps...>;

    using PostProcessingResultVariantT =
        std::invoke_result_t<Postprocessor, std::vector<EvaluationStepResultVariantT>>::value_type;

    using StepResultVariantsT =
        detail::meta::make_step_variants_t<SuccessCollationResultT, EvalSteps...>;

    using StepBuffersT = detail::meta::tuple_of_vector_from_tuple_t<StepResultVariantsT>;

    using LastStepResultVariantT =
        typename std::tuple_element_t<(std::tuple_size_v<StepResultVariantsT> - 1),
                                      StepBuffersT>::value_type;

    static_assert(std::same_as<EvaluationStepResultVariantT, LastStepResultVariantT>);

    using ResultT = EvalResult<PostProcessingResultVariantT, FailedCollationResultT>;

    // Compile time constants
    static constexpr size_t EXPECTED_HIT_GROUP_COUNT = 5;
    static constexpr size_t stepCount{sizeof...(EvalSteps)};
    static constexpr size_t maxStepIndex = stepCount - 1;

    // Dynamic members
    Preprocessor preprocessor;
    Postprocessor postprocessor;
    std::tuple<EvalSteps...> evalSteps;

    StepBuffersT stepBuffers{};
    std::vector<FailedCollationResultT> failureBuffer;

   public:
    constexpr ReadGroupEvaluator(Preprocessor preprocessor, Postprocessor postprocessor,
                                 EvalSteps... steps) noexcept
        : preprocessor(std::move(preprocessor)),
          postprocessor(std::move(postprocessor)),
          evalSteps(std::move(steps)...) {
        std::apply([](auto&... vec) { (vec.reserve(EXPECTED_HIT_GROUP_COUNT), ...); }, stepBuffers);

        failureBuffer.reserve(EXPECTED_HIT_GROUP_COUNT);
    }

    constexpr void clearAll() noexcept {
        std::apply([](auto&... vec) { (vec.clear(), ...); }, stepBuffers);

        failureBuffer.clear();
    }

    [[nodiscard]] auto evaluate(ReadGroup&& readGroup) noexcept -> ResultT {
        clearAll();

        auto& outBuffer = std::get<0>(stepBuffers);

        // 1) Partition items of collation result
        for (auto&& collationResult : preprocessor(std::move(readGroup))) {
            std::visit(
                [&](auto&& constructionResult) {
                    using CollationResultT = std::decay_t<decltype(constructionResult)>;

                    if constexpr (CollationResultT::isFailed()) {
                        failureBuffer.push_back(std::forward<CollationResultT>(constructionResult));
                    } else {
                        outBuffer.emplace_back(std::forward<CollationResultT>(constructionResult));
                    }
                },
                collationResult);
        }

        if constexpr (stepCount > 0) {
            applySteps<0>();
        }

        constexpr std::size_t bufferSize = std::tuple_size_v<decltype(stepBuffers)>;
        auto postprocessedResults = postprocessor(std::move(std::get<bufferSize - 1>(stepBuffers)));

        return {std::move(postprocessedResults), std::move(failureBuffer)};
    }

   private:
    template <std::size_t Index>
    auto applySteps() noexcept {
        auto& step = std::get<Index>(evalSteps);
        auto& inBuffer = std::get<Index>(stepBuffers);
        auto& outBuffer = std::get<Index + 1>(stepBuffers);

        using outBufferValueT =
            std::remove_cvref_t<std::ranges::range_value_t<decltype(outBuffer)>>;

        for (auto&& ctxVariant : inBuffer) {
            std::visit(
                [&](auto&& ctx) {
                    using CtxT = decltype(ctx);
                    if constexpr (!requires { step(std::forward<CtxT>(ctx)); }) {
                        static_assert(detail::meta::is_variant_member_v<std::remove_cvref_t<CtxT>,
                                                                        outBufferValueT> ||
                                          std::is_same_v<void, decltype(step)>,
                                      "Cannot forward unchanged type to out buffer.");
                        outBuffer.emplace_back(std::forward<CtxT>(ctx));
                    } else {
                        processStepResult(step(std::forward<CtxT>(ctx)), outBuffer);
                    }
                },
                std::move(ctxVariant));
        }

        if constexpr (Index < maxStepIndex) {
            applySteps<Index + 1>();
        }
    }

    template <typename R, typename Buffer>
    auto processStepResult(R&& result, Buffer& buffer) noexcept -> void {
        using BufferValueT = std::remove_cvref_t<std::ranges::range_value_t<decltype(buffer)>>;

        if constexpr (is_variant_v<std::remove_cvref_t<R>>) {
            std::visit(
                [&](auto&& inner) {
                    processStepResult(std::forward<decltype(inner)>(inner), buffer);
                },
                std::forward<R>(result));
        } else if constexpr (!std::ranges::range<R>) {
            static_assert(detail::meta::is_variant_member_v<R, BufferValueT>,
                          "Cannot forward step result type to out buffer.");
            buffer.emplace_back(std::forward<R>(result));
        } else {
            static_assert(
                detail::meta::is_variant_member_v<std::ranges::range_value_t<R>, BufferValueT>,
                "Cannot forward step result type to out buffer.");
            std::ranges::move(result, std::back_inserter(buffer));
        }
    }
};

}  // namespace pipelines::detect
