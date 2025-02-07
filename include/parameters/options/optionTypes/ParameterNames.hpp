#pragma once

// Standard
#include <optional>
#include <string>

struct ParameterOptionNames {
    std::optional<char> shortName;
    std::string longName;

    [[nodiscard]] constexpr auto optionsName() const noexcept -> std::string {
        return shortName.has_value() ? longName + "," + std::string(1, shortName.value())
                                     : longName;
    }
};
