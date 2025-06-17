#pragma once

// Standard
#include <cstddef>
#include <optional>
#include <string_view>
#include <tuple>

// Internal
#include "DefaultedParameterOption.hpp"

struct PostprocessOptions {
    static constexpr std::string_view optionsDescription = "Postprocess Pipeline";

    static constexpr DefaultedParameterOption<float, float{0.9}> minSegmentFractionOverlap{
        {.shortName = std::nullopt, .longName = "intfrac"},
        "minimal fractional overlap of each interaction arm between two interactions to be considered the same interaction"sv};

    static constexpr auto allOptions = std::make_tuple(minSegmentFractionOverlap);
};
