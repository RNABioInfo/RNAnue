#pragma once

// Internal
#include "ArithmeticConcept.hpp"
#include "ValueBounds.hpp"

/**
 * @brief Concept to ensure that arithmetic parameters satisfy default value bounds.
 *
 * @tparam T The arithmetic type.
 * @tparam defaultValue The default value.
 * @tparam bounds The valid value bounds.
 */
template <typename T, T defaultValue, ValueBounds<T> bounds>
concept SatisfiesArithmeticParameterConstraints =
    arithmetic<T> && (defaultValue >= bounds.lowerBound) && (defaultValue <= bounds.upperBound);
