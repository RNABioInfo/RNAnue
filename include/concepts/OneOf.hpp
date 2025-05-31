#pragma once

#include <concepts>

template <typename T, typename... Ts>
concept one_of = (std::same_as<T, Ts> or ...);
