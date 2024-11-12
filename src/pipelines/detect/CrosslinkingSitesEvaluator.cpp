#include "CrosslinkingSitesEvaluator.hpp"

// Standard
#include <optional>

// seqan3
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/utility/all.hpp>
#include <utility>

// Internal
#include "Logger.hpp"

namespace pipelines::detect {

auto CrosslinkingSitesEvaluator::evaluate(std::span<const seqan3::dna5> sequence1,
                                          std::span<const seqan3::dna5> sequence2,
                                          const std::vector<seqan3::dot_bracket3> &dotbracket)
    -> std::optional<CrosslinkingSitesEvaluator::Result> {
    if (sequence1.empty() || sequence2.empty() || dotbracket.empty()) [[unlikely]] {
        Logger::log(LogLevel::WARNING, "Empty input sequences or dot-bracket vector!");
        return std::nullopt;
    }

    auto pairedPositions = getPairedNucleotidePositions(dotbracket, sequence1.size());

    if (pairedPositions.empty()) {
        return Result{};
    }

    // Sort pairing sites by first position
    std::ranges::sort(pairedPositions,
                      [](const NucleotidePairPositions &lhs, const NucleotidePairPositions &rhs) {
                          return lhs.first < rhs.first;
                      });

    std::vector<NucleotidePairPositions> intraFirstSequence;
    std::vector<NucleotidePairPositions> intraSecondSequence;
    std::vector<NucleotidePairPositions> interSequences;
    intraFirstSequence.reserve(pairedPositions.size() - 1);
    intraSecondSequence.reserve(pairedPositions.size() - 1);
    interSequences.reserve(pairedPositions.size() - 1);

    // Iterate over all pairs of interaction sites until the second last pair (last pair
    // has no following pair for window)
    for (size_t i = 0; i < pairedPositions.size() - 1; ++i) {
        const NucleotidePositionsWindow window = {pairedPositions[i], pairedPositions[i + 1]};

        auto interactionWindowOpt = getNucleotideWindows(sequence1, sequence2, window);

        if (!interactionWindowOpt) {
            continue;
        }

        InteractionWindow &interactionWindow = interactionWindowOpt.value();

        const NucleotideWindowPair nucleotidePairs = {interactionWindow.forwardWindowNucleotides,
                                                      interactionWindow.reverseWindowNucleotides};

        if (crosslinkingOrientationScheme.contains(nucleotidePairs)) {
            NucleotidePairPositions crosslinkingPositions =
                crosslinkingOrientationScheme.at(nucleotidePairs) ==
                        CrosslinkingOrientation::FORWARD
                    ? std::make_pair(interactionWindow.relativeForwardPositions.first,
                                     interactionWindow.relativeReversePositions.second)
                    : std::make_pair(interactionWindow.relativeForwardPositions.second,
                                     interactionWindow.relativeReversePositions.first);

            if (!interactionWindow.isInterFragment) {
                if (window.second.second < sequence1.size()) {
                    intraFirstSequence.emplace_back(crosslinkingPositions);
                } else {
                    intraSecondSequence.emplace_back(crosslinkingPositions);
                }
            } else {
                interSequences.emplace_back(crosslinkingPositions);
            }
        }
    }

    std::string dotbracketString =
        dotbracket | seqan3::views::to_char | seqan3::ranges::to<std::string>();

    // Dotbracket cannot represent cofolds & fragment intersection -> needs to be reintroduced
    dotbracketString[sequence1.size()] = '&';

    return Result{{intraFirstSequence, intraSecondSequence}, interSequences, dotbracketString};
}

auto CrosslinkingSitesEvaluator::getPairedNucleotidePositions(
    const std::vector<seqan3::dot_bracket3> &dotbracket, const size_t breakpoint)
    -> std::vector<NucleotidePairPositions> {
    std::vector<size_t> openPos;

    std::vector<NucleotidePairPositions> interactionPositions;
    interactionPositions.reserve(dotbracket.size() / 2);

    size_t index = 0;
    bool breakpointReached = false;
    for (const auto &element : dotbracket) {
        if (index == breakpoint && !breakpointReached) {
            breakpointReached = true;
            continue;
        }

        if (element.is_pair_open()) {
            openPos.push_back(index);
        } else if (element.is_pair_close()) {
            const size_t open = openPos.back();
            openPos.pop_back();
            interactionPositions.emplace_back(open, index);
        }
        ++index;
    }

    assert(openPos.empty());

    return interactionPositions;
}

auto CrosslinkingSitesEvaluator::getNucleotideWindows(std::span<const seqan3::dna5> sequence1,
                                                      std::span<const seqan3::dna5> sequence2,
                                                      NucleotidePositionsWindow positionsPair)
    -> std::optional<InteractionWindow> {
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

    // The first pair is in the first sequence and the second pair is in the second sequence
    if (firstPairInSequence1 && !secondPairInSequence1) {
        size_t relativeReverseFirst = reversePair.first - sequence1Length;
        size_t relativeReverseSecond = reversePair.second - sequence1Length;

        return InteractionWindow{
            .forwardWindowNucleotides =
                seqan3::dna5_vector{sequence1[forwardPair.first], sequence1[forwardPair.second]},
            .reverseWindowNucleotides = seqan3::dna5_vector{sequence2[relativeReverseFirst],
                                                            sequence2[relativeReverseSecond]},
            .relativeForwardPositions = forwardPair,
            .relativeReversePositions = std::make_pair(relativeReverseFirst, relativeReverseSecond),
            .isInterFragment = true};
    }

    // Both pairs are in the first sequence
    if (firstPairInSequence1 && secondPairInSequence1) {
        return InteractionWindow{
            .forwardWindowNucleotides =
                seqan3::dna5_vector{sequence1[forwardPair.first], sequence1[forwardPair.second]},
            .reverseWindowNucleotides =
                seqan3::dna5_vector{sequence1[reversePair.first], sequence1[reversePair.second]},
            .relativeForwardPositions = forwardPair,
            .relativeReversePositions = reversePair,
            .isInterFragment = false};
    }

    // Both pairs are in the second sequence
    if (!firstPairInSequence1 && !secondPairInSequence1) {
        size_t relativeForwardFirst = forwardPair.first - sequence1Length;
        size_t relativeForwardSecond = forwardPair.second - sequence1Length;
        size_t relativeReverseFirst = reversePair.first - sequence1Length;
        size_t relativeReverseSecond = reversePair.second - sequence1Length;
        return InteractionWindow{
            .forwardWindowNucleotides = seqan3::dna5_vector{sequence2[relativeForwardFirst],
                                                            sequence2[relativeForwardSecond]},
            .reverseWindowNucleotides = seqan3::dna5_vector{sequence2[relativeReverseFirst],
                                                            sequence2[relativeReverseSecond]},
            .relativeForwardPositions = std::make_pair(relativeForwardFirst, relativeForwardSecond),
            .relativeReversePositions = std::make_pair(relativeReverseFirst, relativeReverseSecond),
            .isInterFragment = false};
    }

    // Due to sorting of base pairs first pair in second sequence and second pair in first
    // sequence is not possible
    return std::nullopt;
}

[[nodiscard]] auto CrosslinkingSitesEvaluator::Result::getTotalCrosslinkingCount() const -> size_t {
    size_t intraCrosslinkingCount = 0;
    for (const auto &intraCrosslinking : intraCrosslinkingSites) {
        intraCrosslinkingCount += intraCrosslinking.size();
    }
    return intraCrosslinkingCount + interCrosslinkingSites.size();
}

[[nodiscard]] auto CrosslinkingSitesEvaluator::Result::getIntraSequenceCrosslinking(
    const size_t fragmentIndex) const -> std::string {
    std::string crosslinkingString;
    for (const auto &crosslinking : intraCrosslinkingSites[fragmentIndex]) {
        crosslinkingString += "[" + std::to_string(crosslinking.first) + "," +
                              std::to_string(crosslinking.second) + "]";
    }
    return crosslinkingString;
}

[[nodiscard]] auto CrosslinkingSitesEvaluator::Result::getInterSequenceCrosslinking(
    const size_t fragmentIndex) const -> std::string {
    std::string crosslinkingString;
    for (const auto &crosslinking : interCrosslinkingSites) {
        crosslinkingString += "[" +
                              (fragmentIndex == 0UL ? std::to_string(crosslinking.first)
                                                    : std::to_string(crosslinking.second)) +
                              "]";
    }
    return crosslinkingString;
}

[[nodiscard]] auto CrosslinkingSitesEvaluator::Result::operator==(
    const CrosslinkingSitesEvaluator::Result &other) const -> bool {
    return this->interCrosslinkingSites == other.interCrosslinkingSites &&
           this->intraCrosslinkingSites == other.intraCrosslinkingSites &&
           this->dotbracket == other.dotbracket;
}

auto operator<<(std::ostream &outputStream, const CrosslinkingSitesEvaluator::Result &result)
    -> std::ostream & {
    outputStream << "Intra-sequence crosslinking: ";
    for (size_t i = 0; i < result.intraCrosslinkingSites.size(); ++i) {
        outputStream << "Fragment " << i << ": " << result.getIntraSequenceCrosslinking(i) << " ";
    }
    outputStream << "\n";
    outputStream << "Inter-sequence crosslinking: ";

    for (size_t i = 0; i < result.interCrosslinkingSites.size(); ++i) {
        outputStream << result.getInterSequenceCrosslinking(i);
    }

    outputStream << "\n";
    outputStream << "Dot-bracket: " << result.dotbracket;
    return outputStream;
}

}  // namespace pipelines::detect
