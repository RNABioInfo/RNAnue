#pragma once

// Standard
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <ostream>
#include <span>
#include <string>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/to_char.hpp>

namespace pipelines::detect {

using namespace seqan3::literals;

class CrosslinkingSitesEvaluator {
   public:
    using NucleotidePairPositions = std::pair<size_t, size_t>;
    using NucleotidePositionsWindow = std::pair<NucleotidePairPositions, NucleotidePairPositions>;

    struct Result;

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

struct CrosslinkingSitesEvaluator::Result {
    Result(const std::vector<std::vector<NucleotidePairPositions>> &intraCrosslinkingSites,
           const std::vector<NucleotidePairPositions> &interCrosslinkingSites,
           std::string dotbracket)
        : intraCrosslinkingSites(intraCrosslinkingSites),
          interCrosslinkingSites(interCrosslinkingSites),
          dotbracket(std::move(dotbracket)) {}

    Result() = default;

    [[nodiscard]] auto getTotalCrosslinkingCount() const -> size_t;

    [[nodiscard]] auto getIntraSequenceCrosslinking(size_t fragmentIndex) const -> std::string;

    [[nodiscard]] auto getInterSequenceCrosslinking(size_t fragmentIndex) const -> std::string;

    [[nodiscard]] auto getDotbracket() const -> const std::string & { return dotbracket; }

    [[nodiscard]] auto operator==(const Result &other) const -> bool;

    friend auto operator<<(std::ostream &outputStream, const Result &result) -> std::ostream &;

   private:
    std::vector<std::vector<NucleotidePairPositions>> intraCrosslinkingSites;
    std::vector<NucleotidePairPositions> interCrosslinkingSites;
    std::string dotbracket;
};

}  // namespace pipelines::detect
