#pragma once

// Standard
#include <cstddef>
#include <numeric>
#include <optional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "CoOptimalPairwiseAligner.hpp"
#include "NucleotidePositions.hpp"
#include "seqan3/alphabet/structure/dot_bracket3.hpp"
#include "seqan3/alphabet/views/to_char.hpp"

namespace dataTypes::HitGroupEvaluation {

struct CrosslinkingResult {
    CrosslinkingResult(
        const std::vector<std::vector<NucleotidePairPositions>> &intraCrosslinkingSites,
        const std::vector<NucleotidePairPositions> &interCrosslinkingSites,
        const std::vector<seqan3::dot_bracket3> &dotbracket, size_t dotbracketStrandEndIndex)
        : intraCrosslinkingSites(intraCrosslinkingSites),
          interCrosslinkingSites(interCrosslinkingSites),
          dotbracket(dotbracket),
          dotbracketStrandEndIndex(dotbracketStrandEndIndex) {}

    CrosslinkingResult() = default;

    [[nodiscard]] constexpr auto getInterCrosslinkingCount() const noexcept -> size_t {
        return interCrosslinkingSites.size();
    };

    [[nodiscard]] constexpr auto getIntraCrosslinkingCount() const noexcept -> size_t {
        return std::accumulate(intraCrosslinkingSites.begin(), intraCrosslinkingSites.end(), 0,
                               [](size_t sum, const auto &crosslinkingSites) {
                                   return sum + crosslinkingSites.size();
                               });
    }

    [[nodiscard]] auto getTotalCrosslinkingCount() const -> size_t {
        size_t intraCrosslinkingCount = 0;
        for (const auto &intraCrosslinking : intraCrosslinkingSites) {
            intraCrosslinkingCount += intraCrosslinking.size();
        }
        return intraCrosslinkingCount + interCrosslinkingSites.size();
    };

    [[nodiscard]] auto getIntraSequenceCrosslinking(size_t fragmentIndex) const -> std::string {
        if (fragmentIndex >= intraCrosslinkingSites.size()) {
            return {"[]"};
        }

        std::string crosslinkingString;
        for (const auto &crosslinking : intraCrosslinkingSites[fragmentIndex]) {
            crosslinkingString += "[" + std::to_string(crosslinking.first) + "," +
                                  std::to_string(crosslinking.second) + "]";
        }
        return crosslinkingString;
    };

    [[nodiscard]] auto getInterSequenceCrosslinking(size_t fragmentIndex) const -> std::string {
        if (fragmentIndex >= intraCrosslinkingSites.size()) {
            return {"[]"};
        }
        std::string crosslinkingString;
        for (const auto &crosslinking : interCrosslinkingSites) {
            crosslinkingString += "[" +
                                  (fragmentIndex == 0UL ? std::to_string(crosslinking.first)
                                                        : std::to_string(crosslinking.second)) +
                                  "]";
        }
        return crosslinkingString;
    };

    [[nodiscard]] auto getDotbracketString() const -> std::string {
        std::string dotbracketString =
            dotbracket | seqan3::views::to_char | seqan3::ranges::to<std::string>();
        dotbracketString[dotbracketStrandEndIndex] = '&';

        // Dotbracket cannot represent cofolds & fragment intersection -> needs to be reintroduced
        dotbracketString[dotbracketStrandEndIndex] = '&';

        return dotbracketString;
    }

    [[nodiscard]] auto operator==(const CrosslinkingResult &other) const -> bool {
        return this->interCrosslinkingSites == other.interCrosslinkingSites &&
               this->intraCrosslinkingSites == other.intraCrosslinkingSites &&
               this->dotbracket == other.dotbracket;
    };

    friend auto operator<<(std::ostream &outputStream, const CrosslinkingResult &result)
        -> std::ostream &;

   private:
    std::vector<std::vector<NucleotidePairPositions>> intraCrosslinkingSites;
    std::vector<NucleotidePairPositions> interCrosslinkingSites;
    std::vector<seqan3::dot_bracket3> dotbracket;
    size_t dotbracketStrandEndIndex;
};

inline auto operator<<(std::ostream &outputStream, const CrosslinkingResult &result)
    -> std::ostream & {
    outputStream << "Intra-sequence crosslinking: ";
    for (size_t i = 0; i < result.intraCrosslinkingSites.size(); ++i) {
        outputStream << "Fragment " << i << ": " << result.getIntraSequenceCrosslinking(i) << " ";
    }
    outputStream << "\n";
    outputStream << "Inter-sequence crosslinking: ";

    outputStream << result.getInterSequenceCrosslinking(0) << ";";
    outputStream << result.getInterSequenceCrosslinking(1);

    outputStream << "\n";
    outputStream << "Dot-bracket: " << result.getDotbracketString();
    return outputStream;
}

struct HybridizationResult {
    double energy;
    std::optional<CrosslinkingResult> crosslinkingResult;
};

using ComplementarityResult = CoOptimalPairwiseAligner::Result;

}  // namespace dataTypes::HitGroupEvaluation
