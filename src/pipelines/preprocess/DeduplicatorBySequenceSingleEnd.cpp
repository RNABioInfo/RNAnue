#include "DeduplicatorBySequenceSingleEnd.hpp"

// Standard
#include <cstddef>
#include <ranges>
#include <unordered_map>
#include <vector>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>

// Internal
#include "DeduplicationOutput.hpp"
#include "FastqRecord.hpp"
#include "HashDNA5Vector.hpp"
#include "Logger.hpp"
#include "SequenceQualityAlgorithms.hpp"  // NOLINT

namespace pipelines::preprocess {

using namespace dataTypes;

auto DeduplicatorBySequenceSingleEnd::deduplicate(const fs::path& recordsPath)
    -> DeduplicationOutputSingle {
    std::unordered_map<std::vector<seqan3::dna5>, DeduplicationRecordSingleEnd, HashDNA5Vector>
        validRecordIDsBySequence;

    size_t duplicateRecords = 0;

    seqan3::sequence_file_input recordInput{recordsPath};

    for (FastqRecord& record : recordInput) {
        const double meanQuality =
            SequenceQualityAlgorithms::meanQualityScore(record.base_qualities());

        if (validRecordIDsBySequence.find(record.sequence()) == validRecordIDsBySequence.end()) {
            validRecordIDsBySequence.emplace(
                record.sequence(),
                DeduplicationRecordSingleEnd{.recordID = record.id(), .meanQuality = meanQuality});
            continue;
        }

        ++duplicateRecords;

        if (meanQuality > validRecordIDsBySequence[record.sequence()].meanQuality) {
            validRecordIDsBySequence.insert_or_assign(
                record.sequence(),
                DeduplicationRecordSingleEnd{.recordID = record.id(), .meanQuality = meanQuality});
        }
    }

    Logger::log("Duplicate records: ", duplicateRecords,
                "; Unique records: ", validRecordIDsBySequence.size());

    auto validRecordIDsView =
        validRecordIDsBySequence | std::views::values |
        std::views::transform([](DeduplicationRecordSingleEnd& record) { return record.recordID; });

    return DeduplicationOutputSingle{
        .validRecordIDs = {validRecordIDsView.begin(), validRecordIDsView.end()}};
}  // namespace DeduplicatorBySequenceSingleEnd::deduplicate(constfs::path&recordsPath)

}  // namespace pipelines::preprocess
