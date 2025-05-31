#pragma once

#include <string_view>

namespace plotting {

struct FileFormat {
    static constexpr std::string_view png = "png";
    static constexpr std::string_view jpeg = "jpeg";
    static constexpr std::string_view eps = "eps";
    static constexpr std::string_view svg = "svg";
};

}  // namespace plotting
