#include "DeduplicatorBySequencePairedEnd.hpp"

// Standard
#include <cstddef>
#include <filesystem>
#include <ranges>
#include <unordered_map>
#include <utility>

// Internal
#include "DeduplicationOutput.hpp"
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

    seqan3::sequence_file_input recForwardIn{recordsFwd};
    seqan3::sequence_file_input recReverseIn{recordsRev};

    for (auto&& [record1, record2] : seqan3::views::zip(recForwardIn, recReverseIn)) {
        const double meanQuality1 =
            SequenceQualityAlgorithms::meanQualityScore(record1.base_qualities());
        const double meanQuality2 =
            SequenceQualityAlgorithms::meanQualityScore(record2.base_qualities());

        const double meanQuality = (meanQuality1 + meanQuality2) / 2;

        auto key = std::make_pair(record1.sequence(), record2.sequence());

        if (validRecordIDsBySequencePair.find(key) == validRecordIDsBySequencePair.end()) {
            validRecordIDsBySequencePair.emplace(
                key,
                DeduplicationRecordPairedEnd{.recordID = record1.id(), .meanQuality = meanQuality});
            continue;
        }

        ++duplicateRecords;

        if (meanQuality > validRecordIDsBySequencePair[key].meanQuality) {
            validRecordIDsBySequencePair.insert_or_assign(
                key,
                DeduplicationRecordPairedEnd{.recordID = record1.id(), .meanQuality = meanQuality});
        }
    }

    auto validRecordIDsView =
        validRecordIDsBySequencePair | std::views::values |
        std::views::transform([](DeduplicationRecordPairedEnd& deduplicatedRecord) {
            return deduplicatedRecord.recordID;
        });

    Logger::log("Duplicate records: ", duplicateRecords,
                "; Unique records: ", validRecordIDsBySequencePair.size());

    return DeduplicationOutputPaired{
        .validRecordIDs = {validRecordIDsView.begin(), validRecordIDsView.end()}};
}

}  // namespace pipelines::preprocess
