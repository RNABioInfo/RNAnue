#include "DeduplicatorBySequenceSingleEnd.hpp"

// Standard
#include <cstddef>

// Internal
#include "DeduplicationOutput.hpp"
#include "Logger.hpp"
#include "SequenceQualityAlgorithms.hpp"  // NOLINT

namespace pipelines::preprocess {

auto DeduplicatorBySequenceSingleEnd::deduplicate(const fs::path& recordsPath)
    -> DeduplicationOutputSingle {
    std::unordered_map<std::string, DeduplicationRecordSingleEnd> recordsMap;

    size_t duplicateRecords = 0;

    seqan3::sequence_file_input recordInput{recordsPath};

    for (FastqRecord& record : recordInput) {
        const double meanQuality =
            SequenceQualityAlgorithms::meanQualityScore(record.base_qualities());

        const auto charView = record.sequence() | std::views::transform([](seqan3::dna5& base) {
                                  return base.to_char();
                              });
        std::string sequence(charView.begin(), charView.end());

        if (recordsMap.find(sequence) == recordsMap.end()) {
            recordsMap[sequence] = {.record = record, .meanQuality = meanQuality};
            continue;
        }

        ++duplicateRecords;

        if (meanQuality > recordsMap[sequence].meanQuality) {
            recordsMap[sequence] = {.record = record, .meanQuality = meanQuality};
        }
    }

    Logger::log("Duplicate records: ", duplicateRecords);
    Logger::log("Unique records: ", recordsMap.size());

    auto recordView =
        recordsMap | std::views::values |
        std::views::transform([](DeduplicationRecordSingleEnd& record) { return record.record; });

    return DeduplicationOutputSingle{.records = {recordView.begin(), recordView.end()}};
}

}  // namespace pipelines::preprocess
