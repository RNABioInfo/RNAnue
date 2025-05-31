#pragma once

#include <type_traits>
#include <variant>

template <typename T>
struct is_variant : std::false_type {};
template <typename... Ts>
struct is_variant<std::variant<Ts...>> : std::true_type {};

template <typename T>
static constexpr bool is_variant_v = is_variant<std::remove_cvref_t<T>>::value;
