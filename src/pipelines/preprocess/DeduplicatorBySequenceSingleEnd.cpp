#include "DeduplicatorBySequenceSingleEnd.hpp"

// Standard
#include <cstddef>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <vector>

// Internal
#include "DeduplicationOutput.hpp"
#include "HashDNA5Vector.hpp"
#include "Logger.hpp"
#include "SequenceQualityAlgorithms.hpp"  // NOLINT

namespace pipelines::preprocess {

auto DeduplicatorBySequenceSingleEnd::deduplicate(const fs::path& recordsPath)
    -> DeduplicationOutputSingle {
    std::unordered_map<std::vector<seqan3::dna5>, DeduplicationRecordSingleEnd, HashDNA5Vector>
        recordsMap;

    size_t duplicateRecords = 0;

    seqan3::sequence_file_input recordInput{recordsPath};

    for (FastqRecord& record : recordInput) {
        const double meanQuality =
            SequenceQualityAlgorithms::meanQualityScore(record.base_qualities());

        if (recordsMap.find(record.sequence()) == recordsMap.end()) {
            recordsMap.emplace(
                record.sequence(),
                DeduplicationRecordSingleEnd{.record = record, .meanQuality = meanQuality});
            continue;
        }

        ++duplicateRecords;

        if (meanQuality > recordsMap[record.sequence()].meanQuality) {
            recordsMap.insert_or_assign(
                record.sequence(),
                DeduplicationRecordSingleEnd{.record = record, .meanQuality = meanQuality});
        }
    }

    Logger::log("Duplicate records: ", duplicateRecords, "; Unique records: ", recordsMap.size());

    auto recordView =
        recordsMap | std::views::values |
        std::views::transform([](DeduplicationRecordSingleEnd& record) { return record.record; });

    return DeduplicationOutputSingle{.records = {recordView.begin(), recordView.end()}};
}

}  // namespace pipelines::preprocess
