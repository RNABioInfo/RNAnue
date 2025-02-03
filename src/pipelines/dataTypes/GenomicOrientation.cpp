#include "GenomicOrientation.hpp"

// Standard
#include <ios>
#include <istream>
#include <ostream>
#include <string>

namespace dataTypes {

auto operator>>(std::istream &input, GenomicOrientation &orientation) -> std::istream & {
    std::string token;
    input >> token;

    if (token == "both") {
        orientation = GenomicOrientation::BOTH;
    } else if (token == "opposite") {
        orientation = GenomicOrientation::OPPOSITE;
    } else if (token == "same") {
        orientation = GenomicOrientation::SAME;
    } else {
        input.setstate(std::ios_base::failbit);
    }
    return input;
}

auto operator<<(std::ostream &output, GenomicOrientation orientation) -> std::ostream & {
    switch (orientation) {
        case GenomicOrientation::BOTH:
            output << "both";
            break;
        case GenomicOrientation::OPPOSITE:
            output << "opposite";
            break;
        case GenomicOrientation::SAME:
            output << "same";
            break;
        default:
            output << "unknown";
            break;
    }

    return output;
}

}  // namespace dataTypes
