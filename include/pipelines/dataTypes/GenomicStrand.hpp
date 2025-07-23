#pragma once

#include <format>
#include <string>
#include <utility>
namespace dataTypes {
enum GenomicStrand : char { FORWARD = '+', REVERSE = '-', NONE = '*' };

inline auto getGenomicStrand(char strand) -> GenomicStrand {
    switch (strand) {
        case '+':
            return FORWARD;
        case '-':
            return REVERSE;
        default:
            return NONE;
    }
}

/**
 * @brief Overload the logical NOT operator for GenomicStrand.
 *
 * @param strand The GenomicStrand to negate.
 * @return GenomicStrand The negated GenomicStrand.
 */
inline auto operator!(GenomicStrand strand) -> GenomicStrand {
    if (strand == FORWARD) {
        return REVERSE;
    }

    if (strand == REVERSE) {
        return FORWARD;
    }

    return NONE;
};

inline static auto toString(const GenomicStrand strand) -> std::string {
    switch (strand) {
        case GenomicStrand::FORWARD:
            return "+";
        case GenomicStrand::REVERSE:
            return "-";
        case GenomicStrand::NONE:
            return ".";
        default:
            std::unreachable();
    }
}

}  // namespace dataTypes

template <>
struct std::formatter<dataTypes::GenomicStrand> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }

    auto format(const dataTypes::GenomicStrand& strand, auto& ctx) const {
        return std::format_to(ctx.out(), "{}", strand);
    }
};
