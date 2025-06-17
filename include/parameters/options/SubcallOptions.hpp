#pragma once

// Standard
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>

// Internal
#include "ParameterOption.hpp"

struct SubcallOptions {
    static constexpr std::string_view optionsDescription = "Subcall Options"sv;

    static constexpr ParameterOption<std::string> subcall{
        {.shortName = std::nullopt, .longName = "subcall"},
        "The subcall to execute. The following subcalls are available: preprocess, align, detect, analyze, postprocess, complete."sv};

    static constexpr auto allOptions = std::make_tuple(subcall);
};
