#include "HybridizationEvaluationStep.hpp"

// Standard
#include <cassert>
#include <cstdlib>
#include <memory>
#include <stdexcept>
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
    const auto sequence2View = record2.sequence() | views::underlying_sequence(record2.flag());

    auto toString = [](const auto &seq) {
        return (seq | seqan3::views::to_char | seqan3::ranges::to<std::string>());
    };

    std::string interactionSeq = toString(sequence1View) + "&" + toString(sequence2View);

    vrna_md_t model;
    vrna_md_set_default(&model);
    model.uniq_ML = 1;  // Required by vrna_subopt's multibranch traceback.
    std::unique_ptr<vrna_fold_compound_t, decltype(&vrna_fold_compound_free)> foldCompound{
        vrna_fold_compound(interactionSeq.c_str(), &model, VRNA_OPTION_DEFAULT | VRNA_OPTION_HYBRID),
        vrna_fold_compound_free};
    if (!foldCompound) throw std::runtime_error("ViennaRNA could not create a fold compound");
    const auto freeSolutions = [](vrna_subopt_sol_s* solutions) {
        if (!solutions) return;
        for (auto* solution = solutions; solution->structure; ++solution) free(solution->structure);
        free(solutions);
    };
    // Keep energy/lexicographic selection of the first zero-band solution.
    std::unique_ptr<vrna_subopt_sol_s, decltype(freeSolutions)> result{
        vrna_subopt(foldCompound.get(), 0, 1, nullptr), freeSolutions};
    if (!result || !result->structure) {
        return {.passed = false, .energy = std::nullopt, .crosslinkingResult = std::nullopt};
    }

    auto secondaryStructure = std::string(result->structure) |
                              seqan3::views::char_to<seqan3::dot_bracket3> |
                              seqan3::ranges::to<std::vector>();

    const auto combinedSequenceLength = sequence1View.size() + sequence2View.size();

    if (secondaryStructure.size() != (combinedSequenceLength + 1)) {
        throw std::runtime_error("ViennaRNA returned a structure with an unexpected length");
    }

    const seqan3::dna5_vector sequence1{sequence1View.begin(), sequence1View.end()};
    const seqan3::dna5_vector sequence2{sequence2View.begin(), sequence2View.end()};

    const auto crosslinkingResult = CrosslinkingSitesEvaluator::evaluate(
        sequence1, sequence2, secondaryStructure, config.includeWobbleBasePairsInCrosslinkingSites);


    return {.passed = isPassingFilters(result->energy),
            .energy = result->energy,
            .crosslinkingResult = crosslinkingResult};
}

auto HybridizationEvaluationStep::isPassingFilters(double energy) const -> bool {
    return energy <= config.mfeThreshold;
}
}  // namespace pipelines::detect
