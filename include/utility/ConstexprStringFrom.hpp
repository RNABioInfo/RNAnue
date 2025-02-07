#pragma once

// Standard
#include <cstdint>

namespace utility {
template <uint8_t... digits>
struct positive_to_chars {
    static const char value[];
};
template <uint8_t... digits>
constexpr char positive_to_chars<digits...>::value[] = {('0' + digits)..., 0};

template <uint8_t... digits>
struct negative_to_chars {
    static const char value[];
};
template <uint8_t... digits>
constexpr char negative_to_chars<digits...>::value[] = {'-', ('0' + digits)..., 0};

template <bool neg, uint8_t... digits>
struct to_chars : positive_to_chars<digits...> {};

template <uint8_t... digits>
struct to_chars<true, digits...> : negative_to_chars<digits...> {};

constexpr uint8_t EXPLODE_CONST = 10;

template <bool neg, uintmax_t rem, uint8_t... digits>
struct explode : explode<neg, rem / EXPLODE_CONST, rem % EXPLODE_CONST, digits...> {};

template <bool neg, uint8_t... digits>
struct explode<neg, 0, digits...> : to_chars<neg, digits...> {};

template <typename T>
constexpr auto cabs(T num) -> uintmax_t {
    return (num < 0) ? -num : num;
}
}  // namespace utility

template <typename Integer, Integer num>
struct string_from : utility::explode<(num < 0), utility::cabs(num)> {};
