#pragma once

namespace pipelines::preprocess {

class Deduplicator {
   public:
    Deduplicator() = delete;

    template <typename Config>
    static auto deduplicate(const Config& config) -> Config::OutputType;
};

}  // namespace pipelines::preprocess
