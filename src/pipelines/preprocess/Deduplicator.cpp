#include "Deduplicator.hpp"

// Internal
#include "DeduplicationConfig.hpp"
#include "DeduplicatorBySequencePairedEnd.hpp"
#include "DeduplicatorBySequenceSingleEnd.hpp"

namespace pipelines::preprocess {

template <>
auto Deduplicator::deduplicate(const DeduplicationBySequenceSingleConfig& config)
    -> DeduplicationBySequenceSingleConfig::OutputType {
    return DeduplicatorBySequenceSingleEnd::deduplicate(config.recordsPath);
}

template <>
auto Deduplicator::deduplicate(const DeduplicationBySequencePairedConfig& config)
    -> DeduplicationBySequencePairedConfig::OutputType {
    return DeduplicatorBySequencePairedEnd::deduplicate(config.recordsPathFwd,
                                                        config.recordsPathRev);
}

}  // namespace pipelines::preprocess
