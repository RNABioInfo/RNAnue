#pragma once

// Standard
#include <variant>

template <class... Args>
struct VariantUnionHelper;

template <class... Args1, class... Args2>
struct VariantUnionHelper<std::variant<Args1...>, std::variant<Args2...>> {
    using type = std::variant<Args1..., Args2...>;
};

template <class Variant1, class Variant2>
using VariantUnion = typename VariantUnionHelper<Variant1, Variant2>::type;
