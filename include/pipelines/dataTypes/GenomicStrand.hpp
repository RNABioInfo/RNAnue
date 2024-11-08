#pragma once

namespace dataTypes {
enum GenomicStrand : char { FORWARD = '+', REVERSE = '-' };

inline auto operator!(GenomicStrand strand) -> dataTypes::GenomicStrand {
    return strand == dataTypes::FORWARD ? dataTypes::REVERSE : dataTypes::FORWARD;
};

}  // namespace dataTypes
