#include "SplitRecordsHybridizationEvaluator.hpp"

// seqan3
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/utility/all.hpp>

// Internal
#include "Logger.hpp"

auto SplitRecordsHybridizationEvaluator::evaluate(
    const SplitRecords &splitRecords,
    const SplitRecordsEvaluationParameters::BaseParameters &parameters)
    -> std::optional<SplitRecordsHybridizationEvaluator::Result> {
    const seqan3::dna5_vector &sequence1 = splitRecords[0].sequence();
    const seqan3::dna5_vector &sequence2 = splitRecords[1].sequence();

    auto toString = [](const auto &seq) {
        return (seq | seqan3::views::to_char | seqan3::ranges::to<std::string>());
    };

    std::string interactionSeq = toString(sequence1) + "&" + toString(sequence2);

    vrna_fold_compound_t *foldCompound = vrna_fold_compound(
        interactionSeq.c_str(), nullptr, VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID);
    std::unique_ptr<char[]> structure(new char[interactionSeq.size() + 1]);  // NOLINT
    float mfe = vrna_cofold(interactionSeq.c_str(), structure.get());

    if (mfe > parameters.mfeThreshold) {
        return std::nullopt;
    }

    auto secondaryStructure = std::string(vrna_cut_point_insert(
                                  structure.get(), static_cast<int>(sequence1.size()) + 1)) |
                              seqan3::views::char_to<seqan3::dot_bracket3> |
                              seqan3::ranges::to<std::vector>();

    std::optional<SplitRecordsHybridizationEvaluator::Result::CrosslinkingResult>
        crosslinkingResult = findCrosslinkingSites(sequence1, sequence2, secondaryStructure);

    vrna_fold_compound_free(foldCompound);

    return SplitRecordsHybridizationEvaluator::Result{.energy = mfe,
                                                      .crosslinkingResult = crosslinkingResult};
}

const std::map<SplitRecordsHybridizationEvaluator::NucleotideWindowPair, size_t>
    SplitRecordsHybridizationEvaluator::crosslinkingScoringScheme = {
        {{"TA"_dna5, "AT"_dna5}, 3},  // Preferred pyrimidine cross-linking
        {{"TG"_dna5, "GT"_dna5}, 3},
        {{"TC"_dna5, "CT"_dna5}, 3},
        {{"AT"_dna5, "TA"_dna5}, 3},
        {{"GT"_dna5, "TG"_dna5}, 3},
        {{"CT"_dna5, "TC"_dna5}, 3},
        {{"CA"_dna5, "AC"_dna5}, 2},  // Non-preferred pyrimidine cross-linking
        {{"CG"_dna5, "GC"_dna5}, 2},
        {{"AC"_dna5, "CA"_dna5}, 2},
        {{"GC"_dna5, "CG"_dna5}, 2},
        {{"TA"_dna5, "GT"_dna5}, 1},  // Wobble base pairs cross-linking
        {{"TG"_dna5, "GT"_dna5}, 1},
        {{"GT"_dna5, "TA"_dna5}, 1},
        {{"GT"_dna5, "TG"_dna5}, 1}};

auto SplitRecordsHybridizationEvaluator::findCrosslinkingSites(
    std::span<const seqan3::dna5> sequence1, std::span<const seqan3::dna5> sequence2,
    std::vector<seqan3::dot_bracket3> &dotbracket)
    -> std::optional<SplitRecordsHybridizationEvaluator::Result::CrosslinkingResult> {
    if (sequence1.empty() || sequence2.empty() || dotbracket.empty()) {
        Logger::log(LogLevel::WARNING, "Empty input sequences or dot-bracket vector!");
        return std::nullopt;
    }

    std::vector<size_t> openPos;

    std::vector<NucleotidePairPositions> interactionPositions;
    interactionPositions.reserve(dotbracket.size() / 2);

    for (size_t i = 0; i < dotbracket.size(); ++i) {
        if (dotbracket[i].is_pair_open()) {
            openPos.push_back(i);
        } else if (dotbracket[i].is_pair_close()) {
            int open = openPos.back();
            openPos.pop_back();
            interactionPositions.emplace_back(open, i);
        }
    }

    // Sort pairing sites by first position
    std::ranges::sort(interactionPositions,
                      [](const NucleotidePairPositions &lhs, const NucleotidePairPositions &rhs) {
                          return lhs.first < rhs.first;
                      });

    if (!openPos.empty()) {
        Logger::log(LogLevel::WARNING, "Unpaired bases in dot-bracket notation! (" +
                                           std::to_string(openPos.size()) + ")");
        return std::nullopt;
    }

    if (interactionPositions.empty()) {
        return Result::CrosslinkingResult{.crosslinkingSites = {},
                                          .normCrosslinkingScore = 0,
                                          .preferredCrosslinkingScore = 0,
                                          .nonPreferredCrosslinkingScore = 0,
                                          .wobbleCrosslinkingScore = 0};
    }

    std::vector<NucleotidePositionsWindow> crosslinkingSites;
    crosslinkingSites.reserve(interactionPositions.size() - 1);

    size_t intraCrosslinkingScore = 0;
    size_t interCrosslinkingScore = 0;

    int preferredCrosslinkingCount = 0;
    int nonPreferredCrosslinkingCount = 0;
    int wobbleCrosslinkingCount = 0;

    // Iterate over all pairs of interaction sites until the second last pair (last pair
    // has no following pair for window)
    for (size_t i = 0; i < interactionPositions.size() - 1; ++i) {
        const NucleotidePositionsWindow window =
            std::make_pair(interactionPositions[i], interactionPositions[i + 1]);

        std::optional<InteractionWindow> interactionWindow =
            getContinuosNucleotideWindows(sequence1, sequence2, window);

        if (interactionWindow) {
            InteractionWindow &v_interactionWindow = interactionWindow.value();

            const NucleotideWindowPair nucleotidePairs =
                std::make_pair(v_interactionWindow.forwardWindowNucleotides,
                               v_interactionWindow.reverseWindowNucleotides);

            const auto iterator = crosslinkingScoringScheme.find(nucleotidePairs);
            if (iterator != crosslinkingScoringScheme.end()) {
                const size_t score = iterator->second;
                if (v_interactionWindow.isInterFragment) {
                    interCrosslinkingScore += score;
                    if (score == 3) {
                        preferredCrosslinkingCount++;
                    } else if (score == 2) {
                        nonPreferredCrosslinkingCount++;
                    } else if (score == 1) {
                        wobbleCrosslinkingCount++;
                    }
                } else {
                    intraCrosslinkingScore += score;
                }
                crosslinkingSites.emplace_back(window);
            }
        }
    }

    size_t minSeqLength = std::min(sequence1.size(), sequence2.size());

    if (minSeqLength == 0) {
        Logger::log(LogLevel::WARNING, "Division by zero!");
        return std::nullopt;
    }

    double normCrosslinkingScore =
        static_cast<double>(interCrosslinkingScore) / static_cast<double>(minSeqLength);
    return Result::CrosslinkingResult{
        .crosslinkingSites = crosslinkingSites,
        .normCrosslinkingScore = normCrosslinkingScore,
        .preferredCrosslinkingScore = preferredCrosslinkingCount,
        .nonPreferredCrosslinkingScore = nonPreferredCrosslinkingCount,
        .wobbleCrosslinkingScore = wobbleCrosslinkingCount};
}

auto SplitRecordsHybridizationEvaluator::getContinuosNucleotideWindows(
    std::span<const seqan3::dna5> sequence1, std::span<const seqan3::dna5> sequence2,
    NucleotidePositionsWindow positionsPair)
    -> std::optional<SplitRecordsHybridizationEvaluator::InteractionWindow> {
    std::pair<uint16_t, uint16_t> forwardPair =
        std::make_pair(positionsPair.first.first, positionsPair.second.first);
    std::pair<uint16_t, uint16_t> reversePair =
        std::make_pair(positionsPair.first.second, positionsPair.second.second);

    // Check that both pairs are from a continuous region
    if (forwardPair.first + 1 != forwardPair.second ||
        reversePair.first - 1 != reversePair.second) {
        return std::nullopt;
    }

    size_t sequence1Length = sequence1.size();

    bool forwardPairSplit =
        forwardPair.first < sequence1Length && forwardPair.second >= sequence1Length;
    bool reversePairSplit =
        reversePair.first >= sequence1Length && reversePair.second < sequence1Length;

    if (forwardPairSplit || reversePairSplit) {
        return std::nullopt;
    }

    bool firstPairInSequence1 = forwardPair.second < sequence1Length;
    bool secondPairInSequence1 = reversePair.first < sequence1Length;

    // Both pairs are in the first sequence
    if (firstPairInSequence1 && secondPairInSequence1) {
        return InteractionWindow{
            .forwardWindowNucleotides =
                seqan3::dna5_vector{sequence1[forwardPair.first], sequence1[forwardPair.second]},
            .reverseWindowNucleotides =
                seqan3::dna5_vector{sequence1[reversePair.first], sequence1[reversePair.second]},
            .forwardWindowPositions = forwardPair,
            .reverseWindowPositions = reversePair,
            .isInterFragment = false};
    }

    // Both pairs are in the second sequence
    if (!firstPairInSequence1 && !secondPairInSequence1) {
        return InteractionWindow{
            .forwardWindowNucleotides =
                seqan3::dna5_vector{sequence2[forwardPair.first - sequence1Length - 1],
                                    sequence2[forwardPair.second - sequence1Length - 1]},
            .reverseWindowNucleotides =
                seqan3::dna5_vector{sequence2[reversePair.first - sequence1Length - 1],
                                    sequence2[reversePair.second - sequence1Length - 1]},
            .forwardWindowPositions = forwardPair,
            .reverseWindowPositions = reversePair,
            .isInterFragment = false};
    }

    // The first pair is in the first sequence and the second pair is in the second sequence
    if (firstPairInSequence1 && !secondPairInSequence1) {
        return InteractionWindow{
            .forwardWindowNucleotides =
                seqan3::dna5_vector{sequence1[forwardPair.first], sequence1[forwardPair.second]},
            .reverseWindowNucleotides =
                seqan3::dna5_vector{sequence2[reversePair.first - sequence1Length - 1],
                                    sequence2[reversePair.second - sequence1Length - 1]},
            .forwardWindowPositions = forwardPair,
            .reverseWindowPositions = reversePair,
            .isInterFragment = true};
    }

    // Due to sorting of base pairs first pair in second sequence and second pair in first
    // sequence is not possible
    return std::nullopt;
}
