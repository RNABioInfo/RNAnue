#pragma once

// Standard
#include <cstdint>
#include <istream>
#include <ostream>

// Internal
#include "GenomicStrandSpecificity.hpp"

namespace dataTypes {

class GenomicOrientation {
   public:
    enum Value : std::uint8_t { SAME, OPPOSITE, BOTH };

    GenomicOrientation() = default;
    constexpr GenomicOrientation(Value orientation) : value(orientation) {};

    constexpr operator Value() const { return value; }

    /**
     * @brief Returns the strand specificity for the current GenomicOrientation.
     *
     * This function is a constexpr method which returns GenomicStrandSpecificity::SPECIFIC
     * if the orientation is either SAME or OPPOSITE, and GenomicStrandSpecificity::UNSPECIFIC
     * otherwise.
     *
     * @return The corresponding GenomicStrandSpecificity.
     */
    [[nodiscard]] constexpr auto strandSpecificity() const noexcept -> GenomicStrandSpecificity {
        switch (value) {
            case SAME:
            case OPPOSITE:
                return GenomicStrandSpecificity::SPECIFIC;
            default:
                return GenomicStrandSpecificity::UNSPECIFIC;
        }
    };

    /**
     * @brief Creates a GenomicOrientation from a given GenomicStrandSpecificity.
     *
     * This function converts a GenomicStrandSpecificity to its corresponding
     * GenomicOrientation. If the specificity is GenomicStrandSpecificity::SPECIFIC,
     * the returned orientation is GenomicOrientation::SAME. Otherwise, it returns
     * GenomicOrientation::BOTH.
     *
     * @param strandSpecificity The GenomicStrandSpecificity to be converted.
     * @return A corresponding GenomicOrientation value.
     */
    [[nodiscard]] static constexpr auto fromStrandSpecificity(
        const GenomicStrandSpecificity strandSpecificity) -> GenomicOrientation {
        return (strandSpecificity == GenomicStrandSpecificity::SPECIFIC) ? GenomicOrientation::SAME
                                                                         : GenomicOrientation::BOTH;
    };

    explicit operator bool() const = delete;

   private:
    Value value;
};

auto operator>>(std::istream& input, GenomicOrientation& orientation) -> std::istream&;
auto operator<<(std::ostream& output, GenomicOrientation orientation) -> std::ostream&;

}  // namespace dataTypes
