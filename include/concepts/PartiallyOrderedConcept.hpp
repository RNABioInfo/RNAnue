#pragma once

#include <concepts>

template <typename T>
concept PartiallyOrdered = std::equality_comparable<T> && requires(const T& lhs, const T& rhs) {
    { lhs < rhs } -> std::convertible_to<bool>;
    { lhs > rhs } -> std::convertible_to<bool>;
};
