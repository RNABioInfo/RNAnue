#pragma once

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <random>
#include <ranges>
#include <vector>

template <typename T>
[[nodiscard]] auto subsample(const std::vector<T>& input, size_t count, std::mt19937& randomDevice)
    -> std::vector<T> {
    if (input.size() <= count) {
        return input;
    }

    std::vector<T> subsampleData;
    subsampleData.reserve(count);
    std::sample(input.begin(), input.end(), std::back_inserter(subsampleData),
                static_cast<int>(count), randomDevice);

    return subsampleData;
}

template <std::ranges::random_access_range R1, std::ranges::random_access_range R2>
auto subsample(const R1& range1, const R2& range2, size_t count) {
    assert(range1.size() == range2.size());
    std::mt19937_64 rng{std::random_device{}()};

    std::vector<std::size_t> res;
    res.reserve(count);
    for (std::size_t index1 = 0, rangeSize = range1.size(); index1 < rangeSize; ++index1) {
        if (res.size() < count) {
            res.push_back(index1);
        } else {
            // choose a random j in [0..i]
            std::uniform_int_distribution<std::size_t> dist(0, index1);
            auto index2 = dist(rng);
            if (index2 < count) {
                res[index2] = index1;  // replace one element
            }
        }
    }

    std::vector<std::ranges::range_value_t<R1>> subRange1;
    std::vector<std::ranges::range_value_t<R2>> subRange2;
    subRange1.reserve(count);
    subRange2.reserve(count);

    for (auto idx : res) {
        subRange1.push_back(range1[idx]);
        subRange2.push_back(range2[idx]);
    }

    return std::make_pair(std::move(subRange1), std::move(subRange2));
}
