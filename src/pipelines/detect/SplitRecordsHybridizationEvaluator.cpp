#include "SplitRecordsHybridizationEvaluator.hpp"

// seqan3
#include <cassert>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/utility/all.hpp>
#include <vector>

// Internal
#include "CrosslinkingSitesEvaluator.hpp"

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
    float mfe = vrna_cofold(interactionSeq.c_str(), structure.get());

    if (mfe > parameters.mfeThreshold) {
        return std::nullopt;
    }

    auto secondaryStructure = std::string(vrna_cut_point_insert(
                                  structure.get(), static_cast<int>(sequence1.size()) + 1)) |
                              seqan3::views::char_to<seqan3::dot_bracket3> |
                              seqan3::ranges::to<std::vector>();

    const auto crosslinkingResult =
        CrosslinkingSitesEvaluator::evaluate(sequence1, sequence2, secondaryStructure);

    vrna_fold_compound_free(foldCompound);

    return SplitRecordsHybridizationEvaluator::Result{.energy = mfe,
                                                      .crosslinkingResult = crosslinkingResult};
}

}  // namespace pipelines::detect
