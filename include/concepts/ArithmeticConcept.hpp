#pragma once

// Standard
#include <concepts>

template <typename T>
concept arithmetic = std::integral<T> or std::floating_point<T>;
