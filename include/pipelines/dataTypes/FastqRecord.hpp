#pragma once

// Standard
#include <vector>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/record.hpp>

namespace dataTypes {
using FastqFieldTypes =
    seqan3::type_list<std::vector<seqan3::dna5>, std::string, std::vector<seqan3::phred42>>;

using FastqFieldIDs = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::qual>;

using FastqRecord = seqan3::sequence_record<FastqFieldTypes, FastqFieldIDs>;
}  // namespace dataTypes
