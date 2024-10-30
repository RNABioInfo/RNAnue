#pragma once

// Standard
#include <filesystem>
#include <utility>
#include <variant>
#include <vector>

// Internal
#include "FastqRecord.hpp"

namespace pipelines::preprocess {

namespace fs = std::filesystem;

using namespace dataTypes;

class Deduplicator {
   public:
    struct DeduplicationRecord {
        FastqRecord record;
        double quality;
    };

    struct BySequenceConfig {
        BySequenceConfig(fs::path inputPath) : recordsPath(std::move(inputPath)) {};

        auto operator()() const -> std::vector<FastqRecord>;

       private:
        fs::path recordsPath;
    };

    using DeduplicationConfigVariant = std::variant<BySequenceConfig>;

    Deduplicator(DeduplicationConfigVariant config) : config(std::move(config)) {};

    auto deduplicate() -> std::vector<FastqRecord>;

   private:
    struct DeduplicationConfig {
        auto operator()(const BySequenceConfig& config) -> std::vector<FastqRecord> {
            return config();
        };
    };

    DeduplicationConfigVariant config;
};

}  // namespace pipelines::preprocess
