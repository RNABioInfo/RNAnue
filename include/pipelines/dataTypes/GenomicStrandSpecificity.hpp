#pragma once

// Standard
#include <ios>
#include <istream>
#include <ostream>
#include <string>

namespace dataTypes {

enum GenomicStrandSpecificity : bool { SPECIFIC = true, UNSPECIFIC = false };

inline auto operator>>(std::istream& input, GenomicStrandSpecificity& strandSpecificity)
    -> std::istream& {
    std::string token;
    input >> token;

    if (token == "specific") {
        strandSpecificity = GenomicStrandSpecificity::SPECIFIC;
    } else if (token == "unspecific") {
        strandSpecificity = GenomicStrandSpecificity::UNSPECIFIC;
    } else {
        input.setstate(std::ios_base::failbit);
    }

    return input;
};

inline auto operator<<(std::ostream& output, GenomicStrandSpecificity strandSpecificity)
    -> std::ostream& {
    switch (strandSpecificity) {
        case GenomicStrandSpecificity::SPECIFIC:
            output << "specific";
            break;
        case GenomicStrandSpecificity::UNSPECIFIC:
            output << "unspecific";
            break;
        default:
            output << "unknown";
            break;
    }

    return output;
};

}  // namespace dataTypes
