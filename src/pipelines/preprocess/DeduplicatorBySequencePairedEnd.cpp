#include "DeduplicatorBySequencePairedEnd.hpp"

// Standard
#include <cstddef>
#include <filesystem>
#include <ranges>
#include <unordered_map>
#include <utility>

// Internal
#include "DeduplicationOutput.hpp"
#include "CheckedFastqReader.hpp"
#include "Logger.hpp"
#include "SequenceQualityAlgorithms.hpp"
#include "seqan3/alphabet/nucleotide/dna5.hpp"
#include "seqan3/io/sequence_file/input.hpp"
#include "utility/PairHash.hpp"

namespace pipelines::preprocess {

auto DeduplicatorBySequencePairedEnd::deduplicate(const fs::path& recordsFwd,
                                                  const fs::path& recordsRev)
    -> DeduplicationOutputPaired {
    std::unordered_map<std::pair<seqan3::dna5_vector, seqan3::dna5_vector>,
                       DeduplicationRecordPairedEnd, PairHash>
        validRecordIDsBySequencePair;

    size_t duplicateRecords = 0;
    size_t nextOrdinal = 0;

    CheckedFastqPairReader input{recordsFwd, recordsRev};
    while (auto pair = input.next()) {
        auto& [record1, record2] = *pair;
        const size_t ordinal = nextOrdinal++;
        const double meanQuality1 =
            SequenceQualityAlgorithms::meanQualityScore(record1.base_qualities());
        const double meanQuality2 =
            SequenceQualityAlgorithms::meanQualityScore(record2.base_qualities());

        const double meanQuality = (meanQuality1 + meanQuality2) / 2;

        auto key = std::make_pair(record1.sequence(), record2.sequence());

        if (!validRecordIDsBySequencePair.contains(key)) {
            validRecordIDsBySequencePair.emplace(
                key,
                DeduplicationRecordPairedEnd{.recordOrdinal = ordinal, .meanQuality = meanQuality});
            continue;
        }

        ++duplicateRecords;

        if (meanQuality > validRecordIDsBySequencePair[key].meanQuality) {
            validRecordIDsBySequencePair.insert_or_assign(
                key,
                DeduplicationRecordPairedEnd{.recordOrdinal = ordinal, .meanQuality = meanQuality});
        }
    }

    auto validRecordIDsView =
        validRecordIDsBySequencePair | std::views::values |
        std::views::transform([](DeduplicationRecordPairedEnd& deduplicatedRecord) {
            return deduplicatedRecord.recordOrdinal;
        });

    Logger::log("Duplicate records: ", duplicateRecords,
                "; Unique records: ", validRecordIDsBySequencePair.size());

    return DeduplicationOutputPaired{
        .validRecordOrdinals = {validRecordIDsView.begin(), validRecordIDsView.end()}};
}

}  // namespace pipelines::preprocess
