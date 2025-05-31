#pragma once

#include <tuple>
#include <utility>

// Generic for_each for a tuple
template <typename Tuple, typename F>
constexpr void tupleForEach(Tuple&& tuple, F&& func) {
    std::apply(
        [&](auto&&... elems) {
            // fold-expression: calls f(elems) for each elem in the pack
            ((void)std::invoke(func, std::forward<decltype(elems)>(elems)), ...);
        },
        std::forward<Tuple>(tuple));
}
