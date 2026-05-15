#pragma once

// Standard
#include <optional>
#include <string>
#include <string_view>

struct ParameterOptionNames {
    std::optional<char> shortName;
    std::string longName;
    std::string_view inverseLongName{};
    std::string_view inverseDescription{};

    [[nodiscard]] constexpr auto optionsName() const noexcept -> std::string {
        return shortName.has_value() ? longName + "," + std::string(1, shortName.value())
                                     : longName;
    }

    [[nodiscard]] constexpr auto hasInverseName() const noexcept -> bool {
        return !inverseLongName.empty();
    }
};
