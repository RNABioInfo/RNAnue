#pragma once

// Standard
#include <filesystem>
#include <string>

// Internal
#include "DeduplicationOutput.hpp"

namespace pipelines::preprocess {

namespace fs = std::filesystem;

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
        std::string recordID;
        double meanQuality;
    };
};

}  // namespace pipelines::preprocess
