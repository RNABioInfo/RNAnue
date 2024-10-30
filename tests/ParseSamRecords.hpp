#pragma once

// Standard
#include <sstream>
#include <vector>

// seqan3
#include <seqan3/io/sam_file/format_sam.hpp>
#include <seqan3/io/sam_file/input.hpp>

// Classes
#include "SamRecord.hpp"

using namespace dataTypes;

inline auto parseSamRecords(const char* samFileRaw) -> std::vector<SamRecord> {
    using sam_file_input_t =
        seqan3::sam_file_input<seqan3::sam_file_input_default_traits<>, SamFieldIDs>;
    sam_file_input_t samFile(std::istringstream{samFileRaw}, seqan3::format_sam{});

    std::vector<SamRecord> samRecordsVector;

    for (auto& record : samFile) {
        samRecordsVector.push_back(record);
    }
    return samRecordsVector;
}
