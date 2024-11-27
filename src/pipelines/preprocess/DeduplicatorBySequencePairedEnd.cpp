#include "DeduplicatorBySequencePairedEnd.hpp"

// standard
#include <string>

// Internal
#include "DeduplicationOutput.hpp"
#include "SequenceQualityAlgorithms.hpp"  // NOLINT
#include "utility/PairHash.hpp"

namespace pipelines::preprocess {

auto DeduplicatorBySequencePairedEnd::deduplicate(const fs::path& recordsFwd,
                                                  const fs::path& recordsRev)
    -> DeduplicationOutputPaired {
    std::unordered_map<std::pair<std::string, std::string>, DeduplicationRecordPairedEnd, PairHash>
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

        const auto charView1 = record1.sequence() | std::views::transform([](seqan3::dna5& base) {
                                   return base.to_char();
                               });
        const auto charView2 = record2.sequence() | std::views::transform([](seqan3::dna5& base) {
                                   return base.to_char();
                               });
        std::string sequence1(charView1.begin(), charView1.end());
        std::string sequence2(charView2.begin(), charView2.end());

        auto key = std::make_pair(sequence1, sequence2);

        if (recordsMap.find(key) == recordsMap.end()) {
            recordsMap[key] = {
                .recordFwd = record1, .recordRev = record2, .meanQuality = meanQuality};
            continue;
        }

        ++duplicateRecords;

        if (meanQuality > recordsMap[key].meanQuality) {
            recordsMap[key] = {
                .recordFwd = record1, .recordRev = record2, .meanQuality = meanQuality};
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

    Logger::log("Duplicate records: ", duplicateRecords);
    Logger::log("Unique records: ", recordsMap.size());

    auto pairView =
        seqan3::views::zip(recordFwdView, recordRevView) |
        std::views::transform([](auto&& pair) { return std::make_pair(pair.first, pair.second); });

    return DeduplicationOutputPaired{.recordPairs = {pairView.begin(), pairView.end()}};
}

}  // namespace pipelines::preprocess
