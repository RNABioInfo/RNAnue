#pragma once

#include <seqan3/alphabet/quality/phred94.hpp>
#include <seqan3/io/sam_file/input.hpp>

namespace utility {
// Transport BAM records without reducing the original quality-score range.
struct AlignmentInputTraits : seqan3::sam_file_input_default_traits<> {
    using quality_alphabet = seqan3::phred94;
};
// Exclude header_ptr so output writers use their own (possibly remapped) header.
using AlignmentFileInput = seqan3::sam_file_input<
    AlignmentInputTraits,
    seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::ref_id,
                   seqan3::field::ref_offset, seqan3::field::cigar, seqan3::field::mapq,
                   seqan3::field::qual, seqan3::field::flag, seqan3::field::mate,
                   seqan3::field::tags>>;
}  // namespace utility
