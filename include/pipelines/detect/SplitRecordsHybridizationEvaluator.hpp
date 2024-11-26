#pragma once

// Standard
#include <optional>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>

// ViennaRNA
extern "C" {
#include <ViennaRNA/subopt.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/strings.h>
}

// Internal
#include "CrosslinkingSitesEvaluator.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"

namespace pipelines::detect {
using namespace seqan3::literals;
using namespace dataTypes;

class SplitRecordsHybridizationEvaluator {
   public:
    struct Result {
        double energy;
        std::optional<pipelines::detect::CrosslinkingSitesEvaluator::Result> crosslinkingResult;
    };

    SplitRecordsHybridizationEvaluator() = delete;
    ~SplitRecordsHybridizationEvaluator() = delete;
    SplitRecordsHybridizationEvaluator(const SplitRecordsHybridizationEvaluator &) = delete;
    auto operator=(const SplitRecordsHybridizationEvaluator &)
        -> SplitRecordsHybridizationEvaluator & = delete;
    SplitRecordsHybridizationEvaluator(SplitRecordsHybridizationEvaluator &&) = delete;
    auto operator=(SplitRecordsHybridizationEvaluator &&)
        -> SplitRecordsHybridizationEvaluator & = delete;

    static auto evaluate(const SplitRecords &splitRecords,
                         const SplitRecordsEvaluationParameters::BaseParameters &parameters)
        -> std::optional<Result>;
};

}  // namespace pipelines::detect
