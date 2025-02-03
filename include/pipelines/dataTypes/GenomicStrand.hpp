#pragma once

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

}  // namespace dataTypes
