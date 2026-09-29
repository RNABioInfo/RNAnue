#pragma once

// Standard
#include <optional>
#include <string>
#include <string_view>

struct ParameterOptionNames {
    std::optional<char> shortName;
    // Names are string literals. A view keeps constexpr options independent of
    // std::string's small-string capacity, including longer snake_case names.
    std::string_view longName;
    std::string_view inverseLongName{};
    std::string_view inverseDescription{};

    [[nodiscard]] constexpr auto optionsName() const noexcept -> std::string {
        return shortName.has_value() ? std::string{longName} + "," + std::string(1, shortName.value())
                                     : std::string{longName};
    }

    [[nodiscard]] constexpr auto hasInverseName() const noexcept -> bool {
        return !inverseLongName.empty();
    }
};
