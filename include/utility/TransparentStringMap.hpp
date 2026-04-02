#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <unordered_map>

namespace util {

struct TransparentStringHash {
    using is_transparent = void;

    [[nodiscard]] auto operator()(std::string_view view) const noexcept -> std::size_t {
        return std::hash<std::string_view>{}(view);
    }

    [[nodiscard]] auto operator()(std::string const& value) const noexcept -> std::size_t {
        return std::hash<std::string_view>{}(std::string_view{value});
    }
};

struct TransparentStringEqual {
    using is_transparent = void;

    [[nodiscard]] auto operator()(std::string const& lhs, std::string const& rhs) const noexcept
        -> bool {
        return lhs == rhs;
    }

    [[nodiscard]] auto operator()(std::string_view lhs, std::string_view rhs) const noexcept
        -> bool {
        return lhs == rhs;
    }

    [[nodiscard]] auto operator()(std::string const& lhs, std::string_view rhs) const noexcept
        -> bool {
        return std::string_view{lhs} == rhs;
    }

    [[nodiscard]] auto operator()(std::string_view lhs, std::string const& rhs) const noexcept
        -> bool {
        return lhs == std::string_view{rhs};
    }
};

using TransparentStringMap =
    std::unordered_map<std::string, std::string, TransparentStringHash, TransparentStringEqual>;

}  // namespace util
