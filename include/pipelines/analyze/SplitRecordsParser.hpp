#pragma once

// Standard
#include <filesystem>
#include <vector>

// seqan3
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/utility/views/chunk.hpp>
// Internal

#include "InteractionCluster.hpp"

namespace pipelines::analyze {

namespace fs = std::filesystem;

struct SplitRecordsParser {
    static auto parse(const fs::path& splitRecordsFilePath) -> std::vector<InteractionCluster> {
        seqan3::sam_file_input splitsIn{splitRecordsFilePath, SamFieldIDs{}};

        std::vector<InteractionCluster> clusters;
        for (auto&& records : splitsIn | seqan3::views::chunk(2)) {
            std::optional<RecordFragment> segment1 =
                RecordFragment::fromSamRecord(*records.begin());
            std::optional<RecordFragment> segment2 =
                RecordFragment::fromSamRecord(*(++records.begin()));

            if (!segment1 || !segment2) [[unlikely]] {
                continue;
            }

            clusters.emplace_back(InteractionCluster::fromRecordFragments(*segment1, *segment2));
        }

        return clusters;
    };
};
}  // namespace pipelines::analyze
