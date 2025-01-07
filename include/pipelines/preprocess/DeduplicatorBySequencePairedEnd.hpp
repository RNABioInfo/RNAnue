#pragma once

// Standard
#include <filesystem>
#include <string>

// Internal
#include "DeduplicationOutput.hpp"

namespace pipelines::preprocess {

namespace fs = std::filesystem;

class DeduplicatorBySequencePairedEnd {
   public:
    DeduplicatorBySequencePairedEnd() = delete;

    static auto hashCombine(size_t lhs, size_t rhs) -> size_t;

    template <class T1, class T2>
    static auto pairHash(const std::pair<T1, T2>& p) -> std::size_t {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);

        return hashCombine(hash1, hash2);
    }

    static auto deduplicate(const fs::path& recordsFwd, const fs::path& recordsRev)
        -> DeduplicationOutputPaired;

   private:
    struct DeduplicationRecordPairedEnd {
        std::string recordID;
        double meanQuality;
    };
};

}  // namespace pipelines::preprocess
