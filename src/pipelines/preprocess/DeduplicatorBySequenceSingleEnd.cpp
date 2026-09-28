#include "DeduplicatorBySequenceSingleEnd.hpp"

// Standard
#include <cstddef>
#include <filesystem>
#include <ranges>
#include <unordered_map>
#include <vector>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>

// Internal
#include "DeduplicationOutput.hpp"
#include "CheckedFastqReader.hpp"
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
    size_t nextOrdinal = 0;

    CheckedFastqReader recordInput{recordsPath};

    while (auto next = recordInput.next()) {
        auto& record = *next;
        const size_t ordinal = nextOrdinal++;
        const double meanQuality =
            SequenceQualityAlgorithms::meanQualityScore(record.base_qualities());

        if (validRecordIDsBySequence.find(record.sequence()) == validRecordIDsBySequence.end()) {
            validRecordIDsBySequence.emplace(
                record.sequence(),
                DeduplicationRecordSingleEnd{.recordOrdinal = ordinal, .meanQuality = meanQuality});
            continue;
        }

        ++duplicateRecords;

        if (meanQuality > validRecordIDsBySequence[record.sequence()].meanQuality) {
            validRecordIDsBySequence.insert_or_assign(
                record.sequence(),
                DeduplicationRecordSingleEnd{.recordOrdinal = ordinal, .meanQuality = meanQuality});
        }
    }

    Logger::log("Duplicate records: ", duplicateRecords,
                "; Unique records: ", validRecordIDsBySequence.size());

    auto validRecordIDsView =
        validRecordIDsBySequence | std::views::values |
        std::views::transform([](DeduplicationRecordSingleEnd& record) { return record.recordOrdinal; });

    return DeduplicationOutputSingle{
        .validRecordOrdinals = {validRecordIDsView.begin(), validRecordIDsView.end()}};
}  // namespace DeduplicatorBySequenceSingleEnd::deduplicate(constfs::path&recordsPath)

}  // namespace pipelines::preprocess
