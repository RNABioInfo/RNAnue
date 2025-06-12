#pragma once

#include <concepts>

template <typename T>
concept NumericType = std::integral<T> || std::floating_point<T>;
