#pragma once

// Standard
#include <concepts>
#include <cstddef>
#include <utility>
#include <vector>

// Internal
#include "EvaluationContext.hpp"
#include "ReadGroup.hpp"
#include "VariantUnion.hpp"

namespace pipelines::detect {

using ConstructedEvaluationContextVariant =
    VariantUnion<EvaluationContextVariant, ErrorEvaluationContextVariant>;

template <typename T>
concept ReadGroupPreprocessor = requires(T preprocessor, ReadGroup&& readGroup) {
    {
        preprocessor(std::move(readGroup))
    } -> std::same_as<std::vector<ConstructedEvaluationContextVariant>>;
};

}  // namespace pipelines::detect
