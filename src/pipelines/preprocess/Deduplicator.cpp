#include "Deduplicator.hpp"

// Standard
#include <array>
#include <functional>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

// seqan3
#include "seqan3/alphabet/nucleotide/dna5.hpp"
#include "seqan3/alphabet/quality/phred42.hpp"
#include "seqan3/io/sequence_file/input.hpp"

// Internal
#include "FastqRecord.hpp"
#include "Utility.hpp"

namespace pipelines::preprocess {

auto Deduplicator::deduplicate() -> std::vector<FastqRecord> {
    return std::visit(DeduplicationConfig(), config);
}

auto Deduplicator::BySequenceConfig::operator()() const -> std::vector<FastqRecord> {
    std::unordered_map<std::string, DeduplicationRecord> recordsMap;

    seqan3::sequence_file_input recordInput{recordsPath};

    for (FastqRecord& record : recordInput) {
        const int sumQualities = fold_left(
            record.base_qualities() |
                std::views::transform([](seqan3::phred42& qual) { return qual.to_phred(); }),
            0, std::plus<>());
        const double meanQuality =
            static_cast<double>(sumQualities) / static_cast<double>(record.base_qualities().size());

        const auto charView = record.sequence() | std::views::transform([](seqan3::dna5& base) {
                                  return base.to_char();
                              });
        std::string sequence(charView.begin(), charView.end());

        if (recordsMap.find(sequence) != recordsMap.end() &&
            recordsMap[sequence].quality >= meanQuality) {
            continue;
        }

        recordsMap[sequence] = {.record = record, .quality = meanQuality};
    }

    auto recordView =
        recordsMap | std::views::values |
        std::views::transform([](DeduplicationRecord& record) { return record.record; });

    return {recordView.begin(), recordView.end()};
}

}  // namespace pipelines::preprocess
