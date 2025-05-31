#pragma once

// Standard
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <span>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/to_char.hpp>

// Internal
#include "HitGroupEvaluationResult.hpp"
#include "NucleotidePositions.hpp"

namespace pipelines::detect {

using namespace seqan3::literals;
using namespace dataTypes;

class CrosslinkingSitesEvaluator {
   public:
    using Result = HitGroupEvaluation::CrosslinkingResult;

    static auto evaluate(std::span<const seqan3::dna5> sequence1,
                         std::span<const seqan3::dna5> sequence2,
                         const std::vector<seqan3::dot_bracket3> &dotbracket,
                         bool includeWobbleBasePairs) -> std::optional<Result>;

   private:
    using NucleotideWindowPair = std::pair<seqan3::dna5_vector, seqan3::dna5_vector>;

    enum class CrosslinkingOrientation : std::uint8_t { FORWARD, REVERSE };

    inline static const std::map<NucleotideWindowPair, CrosslinkingOrientation>
        crosslinkingOrientationSchemeWobble{
            {{"TA"_dna5, "AT"_dna5},
             CrosslinkingOrientation::FORWARD},  // Preferred pyrimidine cross-linking
            {{"AT"_dna5, "TA"_dna5}, CrosslinkingOrientation::REVERSE},
            {{"TA"_dna5, "GT"_dna5},
             CrosslinkingOrientation::FORWARD},  // Wobble base pairs cross-linking
            {{"GT"_dna5, "TA"_dna5}, CrosslinkingOrientation::REVERSE},
            {{"AT"_dna5, "TG"_dna5}, CrosslinkingOrientation::REVERSE},
            {{"TG"_dna5, "AT"_dna5}, CrosslinkingOrientation::FORWARD},
            {{"TG"_dna5, "GT"_dna5}, CrosslinkingOrientation::FORWARD},
            {{"GT"_dna5, "TG"_dna5}, CrosslinkingOrientation::REVERSE}};

    inline static const std::map<NucleotideWindowPair, CrosslinkingOrientation>
        crosslinkingOrientationSchemeNoWobble{
            {{"TA"_dna5, "AT"_dna5},
             CrosslinkingOrientation::FORWARD},  // Preferred pyrimidine cross-linking
            {{"AT"_dna5, "TA"_dna5}, CrosslinkingOrientation::REVERSE}};

    struct InteractionWindow {
        seqan3::dna5_vector forwardWindowNucleotides;
        seqan3::dna5_vector reverseWindowNucleotides;
        NucleotidePairPositions relativeForwardPositions;
        NucleotidePairPositions relativeReversePositions;
        bool isInterFragment;
    };

    static auto getPairedNucleotidePositions(const std::vector<seqan3::dot_bracket3> &dotbracket,
                                             size_t breakpoint)
        -> std::vector<NucleotidePairPositions>;

    static auto getNucleotideWindows(std::span<const seqan3::dna5> sequence1,
                                     std::span<const seqan3::dna5> sequence2,
                                     NucleotidePositionsWindow positionsPair)
        -> std::optional<InteractionWindow>;
};

}  // namespace pipelines::detect
