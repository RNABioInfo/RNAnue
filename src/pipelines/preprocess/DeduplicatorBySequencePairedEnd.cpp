#include "DeduplicatorBySequencePairedEnd.hpp"

// Internal
#include "DeduplicationOutput.hpp"
#include "SequenceQualityAlgorithms.hpp"  // NOLINT
#include "utility/PairHash.hpp"

namespace pipelines::preprocess {

auto DeduplicatorBySequencePairedEnd::deduplicate(const fs::path& recordsFwd,
                                                  const fs::path& recordsRev)
    -> DeduplicationOutputPaired {
    std::unordered_map<std::pair<seqan3::dna5_vector, seqan3::dna5_vector>,
                       DeduplicationRecordPairedEnd, PairHash>
        recordsMap;

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

        if (recordsMap.find(key) == recordsMap.end()) {
            recordsMap.emplace(key, DeduplicationRecordPairedEnd{.recordFwd = std::move(record1),
                                                                 .recordRev = std::move(record2),
                                                                 .meanQuality = meanQuality});
            continue;
        }

        ++duplicateRecords;

        if (meanQuality > recordsMap[key].meanQuality) {
            recordsMap.insert_or_assign(
                key, DeduplicationRecordPairedEnd{.recordFwd = std::move(record1),
                                                  .recordRev = std::move(record2),
                                                  .meanQuality = meanQuality});
        }
    }

    auto recordFwdView = recordsMap | std::views::values |
                         std::views::transform([](DeduplicationRecordPairedEnd& recordPaired) {
                             return recordPaired.recordFwd;
                         });

    auto recordRevView = recordsMap | std::views::values |
                         std::views::transform([](DeduplicationRecordPairedEnd& recordPaired) {
                             return recordPaired.recordRev;
                         });

    assert(std::ranges::distance(recordFwdView) == std::ranges::distance(recordRevView));

    Logger::log("Duplicate records: ", duplicateRecords, "; Unique records: ", recordsMap.size());

    auto pairView =
        seqan3::views::zip(recordFwdView, recordRevView) |
        std::views::transform([](auto&& pair) { return std::make_pair(pair.first, pair.second); });

    return DeduplicationOutputPaired{.recordPairs = {pairView.begin(), pairView.end()}};
}

}  // namespace pipelines::preprocess
