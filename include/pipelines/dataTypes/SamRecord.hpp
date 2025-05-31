#pragma once

// seqan3
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

// seqan3
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/sam_file/record.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

#include "UnderlyingSequence.hpp"

using namespace seqan3::literals;

namespace dataTypes {

using SamFieldTypes =
    seqan3::type_list<std::string, seqan3::sam_flag, std::optional<int32_t>, std::optional<int32_t>,
                      uint8_t, std::vector<seqan3::cigar>, seqan3::dna5_vector,
                      std::vector<seqan3::phred42>, seqan3::sam_tag_dictionary>;

using SamFieldIDs =
    seqan3::fields<seqan3::field::id, seqan3::field::flag, seqan3::field::ref_id,
                   seqan3::field::ref_offset, seqan3::field::mapq, seqan3::field::cigar,
                   seqan3::field::seq, seqan3::field::qual, seqan3::field::tags>;

using SamRecord = seqan3::sam_record<SamFieldTypes, SamFieldIDs>;

/**
 * @brief Get the end position of a SAM record.(half-open)
 * End position of the read exclusive the last postion.
 * @param record The SAM record.
 * @return The end position of the read (zero-based).
 */
inline auto recordEndPosition(const SamRecord& record) -> std::optional<int32_t> {
    const auto start = record.reference_position();

    if (!start.has_value()) {
        return std::nullopt;
    }

    int32_t end = *start;

    for (const auto& cigar : record.cigar_sequence()) {
        if (cigar == 'M'_cigar_operation || cigar == '='_cigar_operation ||
            cigar == 'D'_cigar_operation || cigar == 'N'_cigar_operation ||
            cigar == 'X'_cigar_operation) {
            end += static_cast<int32_t>(get<0>(cigar));
        }
    }

    return end;
}

inline auto softClippedBaseCount(const SamRecord& record) noexcept -> size_t {
    size_t frontSoftClipped = record.cigar_sequence().front() == 'S'_cigar_operation
                                  ? static_cast<size_t>(get<0>(record.cigar_sequence().front()))
                                  : 0;

    size_t backSoftClipped = record.cigar_sequence().back() == 'S'_cigar_operation
                                 ? static_cast<size_t>(get<0>(record.cigar_sequence().front()))
                                 : 0;

    return frontSoftClipped + backSoftClipped;
}

inline auto alignmentLength(const SamRecord& record) noexcept -> size_t {
    return record.sequence().size() - softClippedBaseCount(record);
}

inline auto operator<(const SamRecord& lhs, const SamRecord& rhs) -> bool {
    const auto& lhs_ref_id = lhs.reference_id();
    const auto& rhs_ref_id = rhs.reference_id();

    if (lhs_ref_id != rhs_ref_id) {
        return lhs_ref_id < rhs_ref_id;
    }

    return lhs.reference_position() < rhs.reference_position();
}

inline auto operator>(const SamRecord& lhs, const SamRecord& rhs) -> bool {
    return dataTypes::operator<(rhs, lhs);
}

}  // namespace dataTypes
