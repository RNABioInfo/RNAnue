#pragma once

// Standard
#include <string>

// Internal
#include "ArithmeticConcept.hpp"

template <typename T>
    requires arithmetic<T>
struct ValueBounds {
    T lowerBound;
    T upperBound;

    [[nodiscard]] constexpr auto isValid(const T& value) const noexcept -> bool {
        return value >= lowerBound && value <= upperBound;
    }

    [[nodiscard]] auto description() const noexcept -> std::string {
        return "[" + std::to_string(lowerBound) + "," + std::to_string(upperBound) + "]";
    }
};
