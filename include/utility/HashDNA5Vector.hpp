#pragma once

// Standard
#include <ranges>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>

struct HashDNA5Vector {
    auto operator()(std::vector<seqan3::dna5> const& vec) const -> std::size_t {
        constexpr size_t SHIFT_RIGHT = 16;
        constexpr size_t HASH_CONSTANT_1 = 0x45d9f3b;
        constexpr size_t HASH_CONSTANT_2 = 0x9e3779b9;
        constexpr int SEED_SHIFT_LEFT = 6;
        constexpr int SEED_SHIFT_RIGHT = 2;

        std::size_t seed = vec.size();
        for (auto val :
             vec | std::ranges::views::transform([](auto elem) { return elem.to_rank(); })) {
            val = ((val >> SHIFT_RIGHT) ^ val) * HASH_CONSTANT_1;
            val = ((val >> SHIFT_RIGHT) ^ val) * HASH_CONSTANT_1;
            val = (val >> SHIFT_RIGHT) ^ val;
            seed ^= val + HASH_CONSTANT_2 + (seed << SEED_SHIFT_LEFT) + (seed >> SEED_SHIFT_RIGHT);
        }
        return seed;
    }
};
