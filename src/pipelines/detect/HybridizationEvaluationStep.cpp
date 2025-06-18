#include "HybridizationEvaluationStep.hpp"

// Standard
#include <cassert>
#include <cstdlib>
#include <memory>
#include <optional>
#include <string>
#include <vector>

// viennaRNA
#include <fold_compound.h>
#include <mfe.h>
#include <subopt.h>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/utility/all.hpp>
#include <seqan3/utility/range/to.hpp>

// Internal
#include "CrosslinkingSitesEvaluator.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "SplitRecords.hpp"
#include "seqan3/alphabet/structure/dot_bracket3.hpp"

namespace pipelines::detect {

auto HybridizationEvaluationStep::evaluate(const ChimericRecords &splitRecords) const
    -> HybridizationEvaluationResult {
    const auto &record1 = splitRecords.first();
    const auto &record2 = splitRecords.second();

    const auto sequence1View = record1.sequence() | views::underlying_sequence(record1.flag());
    const auto sequence2View = record2.sequence() | views::underlying_sequence(record1.flag());

    auto toString = [](const auto &seq) {
        return (seq | seqan3::views::to_char | seqan3::ranges::to<std::string>());
    };

    std::string interactionSeq = toString(sequence1View) + "&" + toString(sequence2View);

    vrna_fold_compound_t *foldCompound = vrna_fold_compound(
        interactionSeq.c_str(), nullptr, VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID);

    std::unique_ptr<char[]> structure(new char[interactionSeq.size() + 1]);  // NOLINT

    constexpr int DELTA_MFE = 0;

    std::unique_ptr<vrna_subopt_sol_s, decltype(&free)> result{
        vrna_subopt(foldCompound, DELTA_MFE, 1, nullptr), free};

    if (result == nullptr) {
        vrna_fold_compound_free(foldCompound);
        return {.passed = false, .energy = std::nullopt, .crosslinkingResult = std::nullopt};
    }

    auto secondaryStructure = std::string(result->structure) |
                              seqan3::views::char_to<seqan3::dot_bracket3> |
                              seqan3::ranges::to<std::vector>();

    const auto combinedSequenceLength = sequence1View.size() + sequence2View.size();

    if (secondaryStructure.size() != (combinedSequenceLength + 1)) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
            "Expected size: ", (combinedSequenceLength + 1), ", Got: ", secondaryStructure.size(),
            "\n", secondaryStructure, "\n", std::string(result->structure));
    }

    const seqan3::dna5_vector sequence1{sequence1View.begin(), sequence1View.end()};
    const seqan3::dna5_vector sequence2{sequence2View.begin(), sequence2View.end()};

    const auto crosslinkingResult = CrosslinkingSitesEvaluator::evaluate(
        sequence1, sequence2, secondaryStructure, config.includeWobbleBasePairsInCrosslinkingSites);

    vrna_fold_compound_free(foldCompound);

    return {.passed = isPassingFilters(result->energy),
            .energy = result->energy,
            .crosslinkingResult = crosslinkingResult};
}

auto HybridizationEvaluationStep::isPassingFilters(double energy) const noexcept -> bool {
    return energy <= config.mfeThreshold;
}
}  // namespace pipelines::detect
