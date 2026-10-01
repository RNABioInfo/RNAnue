#pragma once

// Standard
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <ranges>
#include <span>
#include <tuple>

// seqan3
#include <seqan3/alignment/configuration/align_config_band.hpp>
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
   private:
    using SequenceView = std::span<const seqan3::dna5>;

    struct AdapterMatch {
        std::size_t begin;
        std::size_t end;
        std::int64_t score;
    };

    // Preserve the original score threshold, including truncation towards zero.
    static auto adapterMinimumScore(std::size_t overlap, double fraction) -> std::int64_t {
        return static_cast<std::int64_t>(double(overlap) -
                                         ((static_cast<double>(overlap) * fraction) * 2));
    }

    static auto gapsCannotPass(std::size_t adapterLength, std::size_t readLength,
                               std::size_t minOverlap, double fraction) -> bool {
        if (fraction >= 0.5) {
            return false;
        }

        for (std::size_t length = minOverlap; length <= readLength; ++length) {
            const auto threshold = adapterMinimumScore(length, fraction);
            // The threshold is monotone here; longer spans cannot pass either.
            if (threshold > static_cast<std::int64_t>(adapterLength)) {
                break;
            }
            // A gap costs at least 2 + 4 = 6. Relative to the read span a
            // deletion loses >=6 points and an insertion loses >=7 points.
            const auto upperBound = std::min(static_cast<std::int64_t>(adapterLength),
                                             static_cast<std::int64_t>(length) - 6);
            if (upperBound >= threshold) {
                return false;
            }
        }
        return true;
    }

    static auto ungappedScore(SequenceView adapter, SequenceView read, double fraction)
        -> std::optional<std::int64_t> {
        const auto threshold = adapterMinimumScore(read.size(), fraction);
        auto score = static_cast<std::int64_t>(read.size());
        for (std::size_t i = 0; i < read.size(); ++i) {
            // SeqAn's configured simple nucleotide scheme scores equal dna5
            // ranks +1 (including N/N) and unequal ranks -1.
            if (adapter[i] != read[i]) {
                score -= 2;
            }

            if (score < threshold) {
                return std::nullopt;
            }
        }
        return score;
    }

    static void considerMatch(std::optional<AdapterMatch> &best, std::size_t begin, std::size_t end,
                              std::int64_t score) {
        if (!best || score > best->score || (score == best->score && end > best->end)) {
            best = AdapterMatch{.begin = begin, .end = end, .score = score};
        }
    }

    static auto findUngappedAdapter(const Adapter &adapter, SequenceView read,
                                    std::size_t minOverlap) -> std::optional<AdapterMatch> {
        const SequenceView sequence{adapter.sequence};
        const bool fivePrime = adapter.trimmingMode == TrimConfig::Mode::FIVE_PRIME;
        const auto lastStart = read.size() - minOverlap;

        for (std::size_t offset = 0; offset <= lastStart; ++offset) {
            const auto begin = fivePrime ? lastStart - offset : offset;
            const auto remaining = read.size() - begin;
            std::optional<AdapterMatch> best;

            if (remaining >= sequence.size() && sequence.size() >= minOverlap) {
                if (auto score = ungappedScore(sequence, read.subspan(begin, sequence.size()),
                                               adapter.maxMissMatchFraction)) {
                    considerMatch(best, begin, begin + sequence.size(), *score);
                }
            } else if (!fivePrime && remaining < sequence.size()) {
                // A 3' partial adapter must reach the end of the read.
                if (auto score = ungappedScore(sequence.first(remaining), read.subspan(begin),
                                               adapter.maxMissMatchFraction)) {
                    considerMatch(best, begin, read.size(), *score);
                }
            }

            if (fivePrime && begin == 0) {
                // A 5' partial adapter must start at the beginning of the read.
                const auto longest = std::min(sequence.size() - 1, read.size());

                for (std::size_t length = minOverlap; length <= longest; ++length) {
                    if (auto score = ungappedScore(sequence.last(length), read.first(length),
                                                   adapter.maxMissMatchFraction)) {
                        considerMatch(best, 0, length, *score);
                    }
                }
            }

            if (best) {
                return best;
            }
        }
        return std::nullopt;
    }

    static constexpr auto gappedAlignmentConfig =
        seqan3::align_cfg::method_global{} |
        seqan3::align_cfg::scoring_scheme{
            seqan3::nucleotide_scoring_scheme{seqan3::match_score{1}, seqan3::mismatch_score{-1}}} |
        seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-2},
                                           seqan3::align_cfg::extension_score{-4}} |
        seqan3::align_cfg::output_score{};

    class GappedAdapterSearch {
       public:
        GappedAdapterSearch(const Adapter &adapter, SequenceView read, std::size_t minOverlap)
            : adapterSequence{adapter.sequence},
              readSequence{read},
              minOverlap{minOverlap},
              mismatchFraction{adapter.maxMissMatchFraction},
              fivePrime{adapter.trimmingMode == TrimConfig::Mode::FIVE_PRIME},
              maxReadSpan{maximumReadSpan()},
              maxGapBases{scoreDeficit(maxReadSpan) / 4} {}

        auto find() -> std::optional<AdapterMatch> {
            const auto lastStart = readSequence.size() - minOverlap;

            for (std::size_t offset = 0; offset <= lastStart && !finished; ++offset) {
                const auto begin = fivePrime ? lastStart - offset : offset;
                addCandidatesAt(begin);
            }

            if (candidateCount != 0 && !finished) {
                evaluateBatch();
            }

            return best;
        }

       private:
        struct Candidate {
            SequenceView adapter;
            SequenceView read;
            std::size_t begin;
        };

        [[nodiscard]] auto minimumScore(std::size_t readSpan) const -> std::int64_t {
            return adapterMinimumScore(readSpan, mismatchFraction);
        }

        [[nodiscard]] auto scoreDeficit(std::size_t readSpan) const -> std::size_t {
            return static_cast<std::size_t>(static_cast<std::int64_t>(readSpan) -
                                            minimumScore(readSpan));
        }

        [[nodiscard]] auto maximumReadSpan() const -> std::size_t {
            if (mismatchFraction < 0.5) {
                for (std::size_t length = minOverlap; length <= readSequence.size(); ++length) {
                    // Even an all-match alignment cannot score above the adapter length.
                    // The threshold is monotone here, so longer spans cannot pass either.
                    if (minimumScore(length) > static_cast<std::int64_t>(adapterSequence.size())) {
                        return length - 1;
                    }
                }
            }

            return readSequence.size();
        }

        [[nodiscard]] auto passesScoreBounds(const Candidate &candidate) const -> bool {
            const auto adapterLength = candidate.adapter.size();
            const auto readLength = candidate.read.size();
            const auto lengthDifference = adapterLength > readLength ? adapterLength - readLength
                                                                     : readLength - adapterLength;

            if (readLength > maxReadSpan || lengthDifference > maxGapBases) {
                return false;
            }

            auto upperBound = static_cast<std::int64_t>(std::min(adapterLength, readLength));
            if (lengthDifference != 0) {
                // Unequal lengths require at least one gap: open -2, extension -4.
                upperBound -= 2 + 4 * static_cast<std::int64_t>(lengthDifference);
            }

            return upperBound >= minimumScore(readLength);
        }

        static auto seedSearchWindow(SequenceView read, std::size_t seedBegin, std::size_t seedEnd,
                                     std::size_t drift) -> SequenceView {
            const auto windowBegin =
                std::min(read.size(), seedBegin > drift ? seedBegin - drift : 0);
            const auto windowEnd = std::min(read.size(), seedEnd + drift);

            return read.subspan(windowBegin, windowEnd - windowBegin);
        }

        [[nodiscard]] auto hasSurvivingSeed(const Candidate &candidate) const -> bool {
            const auto deficit = scoreDeficit(candidate.read.size());
            const auto maxEdits = deficit / 2;

            // Every edit costs at least two points. With as many possible edits as adapter
            // bases, no nonempty exact seed is guaranteed; skip seed pruning.
            if (maxEdits >= candidate.adapter.size()) {
                return true;
            }

            // In maxEdits + 1 disjoint seeds, at least one must survive an
            // acceptable alignment. Gap bases cost at least four points each.
            const auto seedCount = maxEdits + 1;
            const auto drift = deficit / 4;

            for (std::size_t index = 0; index < seedCount; ++index) {
                const auto seedBegin = index * candidate.adapter.size() / seedCount;
                const auto seedEnd = (index + 1) * candidate.adapter.size() / seedCount;
                const auto seed = candidate.adapter.subspan(seedBegin, seedEnd - seedBegin);
                const auto window = seedSearchWindow(candidate.read, seedBegin, seedEnd, drift);

                if (!std::ranges::search(window, seed).empty()) {
                    return true;
                }
            }

            return false;
        }

        template <typename alignment_config_type>
        void scoreBatch(const alignment_config_type &config) {
            auto pairs = std::span{candidates}.first(candidateCount) |
                         std::views::transform([](const Candidate &candidate) {
                             return std::tuple{candidate.adapter, candidate.read};
                         });
            std::size_t index = 0;

            for (const auto &result : seqan3::align_pairwise(pairs, config)) {
                const auto &candidate = candidates[index++];

                // Starts are visited in trimming order. Finish all candidates
                // at the first passing start before moving to another start.
                if (best && candidate.begin != best->begin) {
                    finished = true;
                    break;
                }

                if (result.score() >= minimumScore(candidate.read.size())) {
                    considerMatch(best, candidate.begin, candidate.begin + candidate.read.size(),
                                  result.score());
                }
            }
        }

        void evaluateBatch() {
            if (maxGapBases <= std::numeric_limits<std::int32_t>::max()) {
                const auto band = static_cast<std::int32_t>(maxGapBases);
                scoreBatch(gappedAlignmentConfig | seqan3::align_cfg::band_fixed_size{
                                                       seqan3::align_cfg::lower_diagonal{-band},
                                                       seqan3::align_cfg::upper_diagonal{band}});
            } else {
                scoreBatch(gappedAlignmentConfig);
            }

            candidateCount = 0;
        }

        void addCandidate(SequenceView adapter, SequenceView read, std::size_t begin) {
            if (best && begin != best->begin) {
                if (candidateCount != 0) {
                    evaluateBatch();
                }

                finished = true;
                return;
            }

            const Candidate candidate{.adapter = adapter, .read = read, .begin = begin};
            if (!passesScoreBounds(candidate) || !hasSurvivingSeed(candidate)) {
                return;
            }

            candidates[candidateCount++] = candidate;
            if (candidateCount == candidates.size()) {
                evaluateBatch();
            }
        }

        void addFullAdapterCandidates(std::size_t begin) {
            const auto shortest = std::max(minOverlap, minimumSpan(adapterSequence.size()));
            const auto longest = std::min(
                {maxReadSpan, adapterSequence.size() + maxGapBases, readSequence.size() - begin});

            for (std::size_t length = shortest; length <= longest && !finished; ++length) {
                addCandidate(adapterSequence, readSequence.subspan(begin, length), begin);
            }
        }

        [[nodiscard]] auto minimumSpan(std::size_t length) const -> std::size_t {
            return length > maxGapBases ? length - maxGapBases : 1;
        }

        void addThreePrimePartialCandidates(std::size_t begin) {
            // A 3' partial adapter must reach the end of the read.
            const auto suffix = readSequence.subspan(begin);
            const auto shortest = minimumSpan(suffix.size());
            const auto longest = std::min(adapterSequence.size() - 1, suffix.size() + maxGapBases);

            for (std::size_t length = shortest; length <= longest && !finished; ++length) {
                addCandidate(adapterSequence.first(length), suffix, begin);
            }
        }

        void addFivePrimePartialCandidates() {
            // A 5' partial adapter must start at the beginning of the read.
            for (std::size_t span = minOverlap; span <= maxReadSpan && !finished; ++span) {
                const auto shortest = minimumSpan(span);
                const auto longest = std::min(adapterSequence.size() - 1, span + maxGapBases);

                for (std::size_t length = shortest; length <= longest && !finished; ++length) {
                    addCandidate(adapterSequence.last(length), readSequence.first(span), 0);
                }
            }
        }

        void addCandidatesAt(std::size_t begin) {
            addFullAdapterCandidates(begin);

            if (finished) {
                return;
            }

            if (!fivePrime) {
                addThreePrimePartialCandidates(begin);
            } else if (begin == 0) {
                addFivePrimePartialCandidates();
            }
        }

        const SequenceView adapterSequence;
        const SequenceView readSequence;
        const std::size_t minOverlap;
        const double mismatchFraction;
        const bool fivePrime;
        const std::size_t maxReadSpan;
        // Deficit = 2*mismatches + 5*insertions + 4*deletions + 2*gap-opens.
        // Therefore an acceptable alignment has at most deficit / 4 gap bases.
        const std::size_t maxGapBases;

        // Bounded batches amortize SeqAn setup without copying sequence data
        // or materializing every possible interval. Results are evaluated lazily.
        std::array<Candidate, 256> candidates{};
        std::size_t candidateCount{0};
        std::optional<AdapterMatch> best;
        bool finished{false};
    };

    [[gnu::noinline]] static auto findGappedAdapter(const Adapter &adapter, SequenceView read,
                                                    std::size_t minOverlap)
        -> std::optional<AdapterMatch> {
        return GappedAdapterSearch{adapter, read, minOverlap}.find();
    }

   public:
    struct TrimWindowedConfig {
        std::size_t windowTrimmingSize;
        std::size_t minMeanWindowPhred;
    };

    RecordTrimmer() = delete;

    /**
     * Remove through the rightmost valid 5' adapter or from the leftmost valid
     * 3' adapter. Partial adapters are accepted only at the corresponding read
     * end
     */
    template <typename record_type>
    static void trimAdapter(const Adapter &adapter, record_type &record,
                            const std::size_t minOverlapTrimming) {
        auto &seq = record.sequence();
        auto &qual = record.base_qualities();
        const auto minOverlap = std::max(std::size_t{1}, minOverlapTrimming);
        if (adapter.sequence.empty() || seq.size() < minOverlap) return;
        const auto match = gapsCannotPass(adapter.sequence.size(), seq.size(), minOverlap,
                                          adapter.maxMissMatchFraction)
                               ? findUngappedAdapter(adapter, SequenceView{seq}, minOverlap)
                               : findGappedAdapter(adapter, SequenceView{seq}, minOverlap);

        if (!match) {
            return;
        }

        if (adapter.trimmingMode == TrimConfig::Mode::FIVE_PRIME) {
            seq.erase(seq.begin(), seq.begin() + match->end);
            qual.erase(qual.begin(), qual.begin() + match->end);
        } else {
            seq.erase(seq.begin() + match->begin, seq.end());
            qual.erase(qual.begin() + match->begin, qual.end());
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
