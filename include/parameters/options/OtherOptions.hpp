#pragma once

// Standard
#include <cstddef>
#include <string>
#include <string_view>
#include <tuple>

// Internal
#include "ParameterOption.hpp"

struct OtherOptions {
    static constexpr std::string_view optionsDescription = "Other"sv;

    static constexpr ParameterOption<bool> printVersion{{.shortName = 'v', .longName = "version"},
                                                        "print version information and exit"sv};

    static constexpr ParameterOption<bool> printHelp{{.shortName = 'h', .longName = "help"},
                                                     "print this help message and exit"sv};

    static constexpr ParameterOption<std::string> configFile{
        {.shortName = 'c', .longName = "config"}, "path to the configuration file"sv};

    static constexpr auto allOptions = std::make_tuple(printVersion, printHelp, configFile);
};
