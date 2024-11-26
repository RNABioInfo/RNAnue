#pragma once

// Standard
#include <filesystem>

// Internal
#include "DeduplicationOutput.hpp"
#include "FastqRecord.hpp"

namespace pipelines::preprocess {

namespace fs = std::filesystem;
using namespace dataTypes;

class DeduplicatorBySequenceSingleEnd {
   public:
    DeduplicatorBySequenceSingleEnd() = delete;
    DeduplicatorBySequenceSingleEnd(const DeduplicatorBySequenceSingleEnd&) = delete;
    auto operator=(const DeduplicatorBySequenceSingleEnd&)
        -> DeduplicatorBySequenceSingleEnd& = delete;
    DeduplicatorBySequenceSingleEnd(DeduplicatorBySequenceSingleEnd&&) = delete;
    auto operator=(DeduplicatorBySequenceSingleEnd&&) -> DeduplicatorBySequenceSingleEnd& = delete;
    ~DeduplicatorBySequenceSingleEnd() = delete;

    static auto deduplicate(const fs::path& recordsPath) -> DeduplicationOutputSingle;

   private:
    struct DeduplicationRecordSingleEnd {
        FastqRecord record;
        double meanQuality;
    };
};

}  // namespace pipelines::preprocess
