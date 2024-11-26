#pragma once

// Standard
#include <cstddef>
#include <utility>

struct PairHash {
    static auto hashCombine(size_t lhs, size_t rhs) -> size_t {
        constexpr size_t SIZE_THRESHOLD = 8;
        constexpr size_t HASH_CONSTANT_64 = 0x517cc1b727220a95;
        constexpr size_t HASH_CONSTANT_32 = 0x9e3779b9;
        constexpr int SHIFT_LEFT = 6;
        constexpr int SHIFT_RIGHT = 2;

        if constexpr (sizeof(size_t) >= SIZE_THRESHOLD) {
            lhs ^= rhs + HASH_CONSTANT_64 + (lhs << SHIFT_LEFT) + (lhs >> SHIFT_RIGHT);
        } else {
            lhs ^= rhs + HASH_CONSTANT_32 + (lhs << SHIFT_LEFT) + (lhs >> SHIFT_RIGHT);
        }
        return lhs;
    }

    template <class T1, class T2>
    auto operator()(const std::pair<T1, T2> &pair) const -> std::size_t {
        auto hash1 = std::hash<T1>{}(pair.first);
        auto hash2 = std::hash<T2>{}(pair.second);

        return hashCombine(hash1, hash2);
    }
};
