#pragma once

// Standard
#include <concepts>

namespace pipelines::detect {

template <typename P>
concept ReadGroupEvaluationPolicy = requires(P const& policy) {
    { policy.makeEvaluator() };

    { policy.minComplementarity() } -> std::convertible_to<double>;
    { policy.splicingTolerance() } -> std::convertible_to<float>;
    { policy.allow_alternative_splicing() } -> std::convertible_to<bool>;
};

}  // namespace pipelines::detect
