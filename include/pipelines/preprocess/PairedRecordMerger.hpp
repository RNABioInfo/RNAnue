#pragma once

// Standard
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <optional>
#include <ranges>
#include <tuple>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/views/complement.hpp>

// Internal
#include "FastqRecord.hpp"

namespace pipelines::preprocess {

using namespace dataTypes;

struct PairedRecordMerger {
    /**
     * @brief Merges two records if they have a sufficient overlap.
     *
     * This function takes two records and checks if they have a sufficient overlap based on the
     * specified minimum overlap value. If the overlap is sufficient, the records are merged into a
     * new record and returned as an optional. Otherwise, std::nullopt is returned.
     *
     * @tparam record_type The type of the records.
     * @param record1 The first record to be merged.
     * @param record2 The second record to be merged.
     * @param minOverlap The minimum required overlap between the records.
     * @return An optional containing the merged record if the overlap is sufficient, otherwise
     * std::nullopt.
     */
    static auto mergeRecordPair(const PairedFastqRecords &records,
                                size_t minOverlapMerge,  // NOLINT
                                double maxMissMatchRateMerge)
        -> std::optional<FastqRecord> {  // NOLINT
        const seqan3::align_cfg::method_global endGapConfig{
            seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
            seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
            seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
            seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}};

        const seqan3::align_cfg::scoring_scheme scoringSchemeConfig{
            seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}}};

        const seqan3::align_cfg::gap_cost_affine gapSchemeConfig{
            seqan3::align_cfg::open_score{-2}, seqan3::align_cfg::extension_score{-4}};

        const auto outputConfig =
            seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_begin_position{} |
            seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_alignment{};

        const auto alignmentConfig =
            endGapConfig | scoringSchemeConfig | gapSchemeConfig | outputConfig;

        const auto &seq1 = records.first.sequence();
        const auto &seq2ReverseComplement =
            records.second.sequence() | std::views::reverse | seqan3::views::complement;

        std::optional<FastqRecord> mergedRecord{std::nullopt};

        for (auto const &result :
             seqan3::align_pairwise(std::tie(seq1, seq2ReverseComplement), alignmentConfig)) {
            const int overlap = static_cast<int>(result.sequence1_end_position()) -
                                static_cast<int>(result.sequence1_begin_position());
            if (overlap < int(minOverlapMerge)) {
                continue;
            }

            const int minScore =
                static_cast<int>(overlap - ((overlap * maxMissMatchRateMerge) * 2));
            if (result.score() >= minScore) {
                mergedRecord = constructMergedRecord(records, result);
            }
        }

        return mergedRecord;
    }

    /**
     * Constructs a merged record by combining two input records with a specified overlap.
     *
     * @tparam record_type The type of the input records.
     * @param record1 The first input record.
     * @param record2 The second input record.
     * @param overlap The length of the overlap between the two records.
     * @return The merged record.
     */
    template <typename result_type>
    static auto constructMergedRecord(const PairedFastqRecords &records,
                                      const seqan3::alignment_result<result_type> &alignmentResult)
        -> FastqRecord {
        const auto &record1Qualities = records.first.base_qualities();
        const auto &record2Qualities = records.second.base_qualities();

        seqan3::dna5_vector mergedSequence{};
        std::vector<seqan3::phred42> mergedQualities{};

        mergedSequence.reserve(records.first.sequence().size() + records.second.sequence().size());
        mergedQualities.reserve(record1Qualities.size() + record2Qualities.size());

        const auto &record2RevCompSequence =
            records.second.sequence() | seqan3::views::complement | std::views::reverse;
        const auto &record2RevCompQualities = record2Qualities | std::views::reverse;

        // 5' overhang of read one that is not in overlap region
        mergedSequence.insert(
            mergedSequence.end(), records.first.sequence().begin(),
            records.first.sequence().begin() + alignmentResult.sequence1_begin_position());
        mergedQualities.insert(
            mergedQualities.end(), records.first.base_qualities().begin(),
            records.first.base_qualities().begin() + alignmentResult.sequence1_begin_position());

        size_t posRecord1 = alignmentResult.sequence1_begin_position();
        size_t posRecord2 = alignmentResult.sequence2_begin_position();

        const auto &[alignmentSeq1, alignmentSeq2] = alignmentResult.alignment();

        for (const auto &[el1, el2] : seqan3::views::zip(alignmentSeq1, alignmentSeq2)) {
            const seqan3::phred42 qual1 = record1Qualities[posRecord1];
            const seqan3::phred42 qual2 = record2RevCompQualities[static_cast<long>(posRecord2)];

            if (el1 == el2) {
                const auto base1 = el1.template convert_to<seqan3::dna5>();
                mergedSequence.insert(mergedSequence.end(), base1);
                mergedQualities.insert(mergedQualities.end(), std::max(qual1, qual2));

                posRecord1++;
                posRecord2++;
                continue;
            }

            if (el1 != seqan3::gap{} && el2 != seqan3::gap{}) {
                if (qual1 < qual2) {
                    const auto base2 = el2.template convert_to<seqan3::dna5>();
                    mergedSequence.insert(mergedSequence.end(), base2);
                    mergedQualities.insert(mergedQualities.end(), qual2);
                } else {
                    const auto base1 = el1.template convert_to<seqan3::dna5>();
                    mergedSequence.insert(mergedSequence.end(), base1);
                    mergedQualities.insert(mergedQualities.end(), qual1);
                }

                posRecord1++;
                posRecord2++;
                continue;
            }

            if (el1 == seqan3::gap{}) {
                const auto base2 = el2.template convert_to<seqan3::dna5>();
                mergedSequence.insert(mergedSequence.end(), base2);
                mergedQualities.insert(mergedQualities.end(), qual2);

                posRecord2++;
                continue;
            }

            if (el2 == seqan3::gap{}) {
                const auto base1 = el1.template convert_to<seqan3::dna5>();
                mergedSequence.insert(mergedSequence.end(), base1);
                mergedQualities.insert(mergedQualities.end(), qual1);

                posRecord1++;
                continue;
            }
        }

        // 5' overhang of read two that is not in overlap region
        mergedSequence.insert(
            mergedSequence.end(),
            record2RevCompSequence.begin() + alignmentResult.sequence2_end_position(),
            record2RevCompSequence.end());
        mergedQualities.insert(
            mergedQualities.end(),
            record2RevCompQualities.begin() + alignmentResult.sequence2_end_position(),
            record2RevCompQualities.end());

        assert(mergedSequence.size() == mergedQualities.size());

        return FastqRecord{std::move(mergedSequence), records.first.id(),
                           std::move(mergedQualities)};
    }
};

}  // namespace pipelines::preprocess
