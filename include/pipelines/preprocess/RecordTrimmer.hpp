#pragma once

// Standard
#include <cstddef>

// seqan3
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/utility/views/slice.hpp>

// Internal
#include "Adapter.hpp"
#include "SequenceQualityAlgorithms.hpp"
#include "TrimConfig.hpp"

using seqan3::operator""_dna5;

namespace pipelines::preprocess {

struct RecordTrimmer {
    struct TrimWindowedConfig {
        std::size_t windowTrimmingSize;
        std::size_t minMeanWindowPhred;
    };

    RecordTrimmer() = delete;

    /**
     * Trims adapter sequences from the given record.
     *
     * @param adapterSequence The adapter sequence to be trimmed.
     * @param record The record from which the adapter sequence will be trimmed.
     * @param trimmingMode The trimming mode to be applied.
     */
    template <typename record_type>
    static void trimAdapter(const Adapter &adapter, record_type &record,
                            const std::size_t minOverlapTrimming) {
        const seqan3::align_cfg::scoring_scheme scoringSchemeConfig{
            seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}}};
        const seqan3::align_cfg::gap_cost_affine gapSchemeConfig{
            seqan3::align_cfg::open_score{-2}, seqan3::align_cfg::extension_score{-4}};
        const auto outputConfig = seqan3::align_cfg::output_score{} |
                                  seqan3::align_cfg::output_end_position{} |
                                  seqan3::align_cfg::output_begin_position{};

        const auto alignment_config = TrimConfig::alignmentConfigFor(adapter.trimmingMode) |
                                      scoringSchemeConfig | gapSchemeConfig | outputConfig;

        const auto seqCopy = record.sequence();
        auto &seq = record.sequence();
        auto &qual = record.base_qualities();

        for (auto const &result :
             seqan3::align_pairwise(std::tie(adapter.sequence, seqCopy), alignment_config)) {
            const int overlap = result.sequence2_end_position() - result.sequence2_begin_position();

            if (overlap < int(minOverlapTrimming)) {
                continue;
            }

            const int minScore = overlap - ((overlap * adapter.maxMissMatchFraction) * 2);

            if (result.score() >= minScore) {
                if (adapter.trimmingMode == TrimConfig::Mode::FIVE_PRIME) {
                    seq.erase(seq.begin(), seq.begin() + result.sequence2_end_position());
                    qual.erase(qual.begin(), qual.begin() + result.sequence2_end_position());
                } else {
                    seq.erase(seq.begin() + result.sequence2_begin_position(), seq.end());
                    qual.erase(qual.begin() + result.sequence2_begin_position(), qual.end());
                }

                break;
            }
        }
    }

    /**
     * @brief Trims trailing polyG bases from the 3' end of a sequence record.
     *
     * This function removes any consecutive 'G' bases from the end of the sequence
     * if they meet the quality threshold. The mean quality threshold is set to a rank of 20
     * using the phred42 scale. If the number of consecutive polyG bases is greater than
     * or equal to 5, they are removed from both the sequence and the base qualities.
     *
     * @tparam record_type The type of the sequence record.
     * @param record The sequence record to trim.
     */
    template <typename record_type>
    static void trim3PolyG(record_type &record, const size_t minPolyGCount) {
        auto &seq = record.sequence();
        auto &qual = record.base_qualities();

        constexpr int qualityThresholdRank = 20;

        auto seqIt = seq.rbegin();
        auto qualIt = qual.rbegin();

        std::size_t sumPhred = 0;
        std::size_t count = 0;

        while (seqIt != seq.rend() && *seqIt == 'G'_dna5) {
            const auto currentPhred = seqan3::to_phred(*qualIt);
            const auto newCount = count + 1;
            const auto average =
                static_cast<double>(sumPhred + currentPhred) / static_cast<double>(newCount);

            if (average >= qualityThresholdRank) {
                sumPhred += currentPhred;
                ++count;
                ++seqIt;
                ++qualIt;
            } else {
                break;
            }
        }

        if (count >= minPolyGCount) {
            seq.erase(seq.end() - count, seq.end());
            qual.erase(qual.end() - count, qual.end());
        }
    }

    /**
     * Trims the windowed quality of a given record.
     *
     * @tparam record_type The type of the record.
     * @param record The record to trim.
     */
    template <typename record_type>
    static void trimWindowedQuality(record_type &record, const TrimWindowedConfig &config) {
        auto trimmingEnd = static_cast<std::ptrdiff_t>(record.sequence().size());

        while ((trimmingEnd - static_cast<std::ptrdiff_t>(config.windowTrimmingSize)) >=
               static_cast<std::ptrdiff_t>(config.windowTrimmingSize)) {
            auto windowQual =
                record.base_qualities() |
                seqan3::views::slice(
                    trimmingEnd - static_cast<std::ptrdiff_t>(config.windowTrimmingSize),
                    trimmingEnd);

            const double meanQuality = SequenceQualityAlgorithms::meanQualityScore(windowQual);

            if (meanQuality >= static_cast<double>(config.minMeanWindowPhred)) {
                break;
            }
            trimmingEnd--;
        }

        record.sequence().erase(record.sequence().begin() + trimmingEnd, record.sequence().end());
        record.base_qualities().erase(record.base_qualities().begin() + trimmingEnd,
                                      record.base_qualities().end());
    }
};

}  // namespace pipelines::preprocess
