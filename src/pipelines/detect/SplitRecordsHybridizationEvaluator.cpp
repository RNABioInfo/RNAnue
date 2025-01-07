#include "SplitRecordsHybridizationEvaluator.hpp"

// Standard
#include <cassert>
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
#include "Logger.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"

namespace pipelines::detect {

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

    constexpr int DELTA_MFE = 0;
    std::unique_ptr<vrna_subopt_sol_s, decltype(&free)> result{
        vrna_subopt(foldCompound, DELTA_MFE, 1, nullptr), free};

    if (result == nullptr || result->energy > parameters.mfeThreshold) {
        vrna_fold_compound_free(foldCompound);
        return std::nullopt;
    }

    auto secondaryStructure = std::string(result->structure) |
                              seqan3::views::char_to<seqan3::dot_bracket3> |
                              seqan3::ranges::to<std::vector>();

    if (secondaryStructure.size() != (sequence1.size() + sequence2.size() + 1)) {
        Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
            "Expected size: ", (sequence1.size() + sequence2.size() + 1),
            ", Got: ", secondaryStructure.size(), "\n", secondaryStructure, "\n",
            std::string(result->structure));
    }

    const auto crosslinkingResult =
        CrosslinkingSitesEvaluator::evaluate(sequence1, sequence2, secondaryStructure,
                                             parameters.includeWobbleBasePairsInCrosslinkingSites);

    vrna_fold_compound_free(foldCompound);

    return SplitRecordsHybridizationEvaluator::Result{.energy = result->energy,
                                                      .crosslinkingResult = crosslinkingResult};
}

}  // namespace pipelines::detect
