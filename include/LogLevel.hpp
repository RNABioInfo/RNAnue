#pragma once

#include <cstdint>
#include <iostream>
#include <istream>
#include <ostream>
#include <string>
#include <type_traits>

enum LogLevel : std::uint8_t { DEBUG, INFO, WARNING, ERROR };

[[nodiscard]] constexpr auto operator==(LogLevel lhs, LogLevel rhs) noexcept -> bool {
    return static_cast<std::underlying_type_t<LogLevel>>(lhs) ==
           static_cast<std::underlying_type_t<LogLevel>>(rhs);
}

[[nodiscard]] constexpr auto operator!=(LogLevel lhs, LogLevel rhs) noexcept -> bool {
    return !(lhs == rhs);
}

[[nodiscard]] constexpr auto operator<(LogLevel lhs, LogLevel rhs) noexcept -> bool {
    return static_cast<std::underlying_type_t<LogLevel>>(lhs) <
           static_cast<std::underlying_type_t<LogLevel>>(rhs);
}

[[nodiscard]] constexpr auto operator>(LogLevel lhs, LogLevel rhs) noexcept -> bool {
    return rhs < lhs;
}

[[nodiscard]] constexpr auto operator<=(LogLevel lhs, LogLevel rhs) noexcept -> bool {
    return !(rhs < lhs);
}

[[nodiscard]] constexpr auto operator>=(LogLevel lhs, LogLevel rhs) noexcept -> bool {
    return !(lhs < rhs);
}

//----------------------------------------------------------------------------
// Unary minus operator for LogLevel
//----------------------------------------------------------------------------
// Boost's implementation (or related SFINAE checks) may try to apply a unary minus
// to the default value. Since LogLevel is an enum class (with no built-in unary-),
// we define one. Note: this conversion is only meant to satisfy Boost’s templates;
// the resulting value is not used in any arithmetic sense.
[[nodiscard]] constexpr auto operator-(LogLevel lvl) noexcept -> LogLevel {
    return static_cast<LogLevel>(-static_cast<std::underlying_type_t<LogLevel>>(lvl));
}

inline auto operator>>(std::istream& input, LogLevel& level) -> std::istream& {
    std::string token;
    input >> token;
    if (token == "debug") {
        level = LogLevel::DEBUG;
    } else if (token == "info") {
        level = LogLevel::INFO;
    } else if (token == "warning") {
        level = LogLevel::WARNING;
    } else if (token == "error") {
        level = LogLevel::ERROR;
    } else {
        input.setstate(std::ios_base::failbit);
    }
    return input;
}

inline auto operator<<(std::ostream& output, LogLevel level) -> std::ostream& {
    switch (level) {
        case LogLevel::DEBUG:
            output << "debug";
            break;
        case LogLevel::INFO:
            output << "info";
            break;
        case LogLevel::WARNING:
            output << "warning";
            break;
        case LogLevel::ERROR:
            output << "error";
            break;
        default:
            output << "unknown";
            break;
    }
    return output;
}
