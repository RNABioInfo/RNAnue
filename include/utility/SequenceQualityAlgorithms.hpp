#pragma once

// Standard
#include <array>
#include <cmath>
#include <concepts>

// seqan3
#include <seqan3/alphabet/quality/phred42.hpp>

// Internal
#include "Utility.hpp"

namespace SequenceQualityAlgorithms {

template <typename T>
concept IterableSequenceQuality =
    std::same_as<typename T::value_type, seqan3::phred42> && requires(T container) {
        container.begin();
        container.end();
    };

using namespace seqan3::literals;

static constexpr double PHRED_SCALE_BASE = 10;
static constexpr const unsigned char ALPHABET_SIZE = seqan3::phred42::alphabet_size;
using PhredToProbabilityTable = std::array<double, ALPHABET_SIZE>;

// NOLINTBEGIN
consteval auto getPhredToProbabilityTable() -> PhredToProbabilityTable {
    PhredToProbabilityTable table{};

    for (unsigned char phred = 0; phred < seqan3::phred42::alphabet_size; ++phred) {
        table[phred] = pow(PHRED_SCALE_BASE, -static_cast<double>(phred) / PHRED_SCALE_BASE);
    }
    return table;
};

static_assert(getPhredToProbabilityTable()[0] == 1.0);
static_assert(getPhredToProbabilityTable()[10] == 0.1);
static_assert(getPhredToProbabilityTable()[20] == 0.01);
static_assert(getPhredToProbabilityTable()[30] == 0.001);
static_assert(getPhredToProbabilityTable()[40] == 0.0001);

static constexpr PhredToProbabilityTable phredToProbabilityTable{getPhredToProbabilityTable()};

constexpr auto probability(const seqan3::phred42 phred) -> double {
    return phredToProbabilityTable[phred.to_phred() + seqan3::phred42::offset_phred];
};

static_assert(probability('!'_phred42) == 1.0);
static_assert(probability('5'_phred42) == 0.01);

template <IterableSequenceQuality T>
constexpr auto meanQualityScore(const T &qualities) -> double {
    typename T::const_iterator iter;

    double qualitySum = 0.0;

    for (iter = qualities.begin(); iter != qualities.end(); ++iter) {
        qualitySum += probability(*iter);
    }

    const double meanErrorProbability =
        static_cast<double>(qualitySum) / static_cast<double>(qualities.size());

    return (-PHRED_SCALE_BASE * std::log10(meanErrorProbability));
};

// std::span does not support const_iterator as of c++20 (remove with c++23)
template <IterableSequenceQuality T>
constexpr auto meanQualityScore(T &qualities) -> double {
    double qualitySum = 0.0;

    for (const auto qual : qualities) {
        qualitySum += probability(qual);
    }

    const double meanErrorProbability = qualitySum / static_cast<double>(qualities.size());

    return (-PHRED_SCALE_BASE * std::log10(meanErrorProbability));
};

static_assert(meanQualityScore(std::array<seqan3::phred42, 1>{'!'_phred42}) == 0.0);
static_assert(meanQualityScore(std::array<seqan3::phred42, 1>{'5'_phred42}) == 20.0);
static_assert(helper::isApproxEqual(
    meanQualityScore(std::array<seqan3::phred42, 2>{'!'_phred42, '5'_phred42}), 2.967086219));

// NOLINTEND
};  // namespace SequenceQualityAlgorithms
