#pragma once

// Standard
#include <deque>
#include <iostream>
#include <variant>

// Internal
#include "SplitRecordsComplementarityEvaluator.hpp"
#include "SplitRecordsEvaluationParameters.hpp"
#include "SplitRecordsHybridizationEvaluator.hpp"

using namespace seqan3::literals;

namespace pipelines::detect {

class SplitRecordsEvaluator {
   public:
    struct EvaluatedSplitRecords {
        SplitRecords splitRecords;
        SplitRecordsComplementarityEvaluator::Result complementarityResult;
        SplitRecordsHybridizationEvaluator::Result hybridizationResult;

        auto operator<(const EvaluatedSplitRecords &other) const -> bool;
        auto operator>(const EvaluatedSplitRecords &other) const -> bool;
    };

    enum class FilterReason { NO_SPLIT_READ, UNMAPPED, SPLICING, COMPLEMENTARITY, HYBRIDIZATION };

    using Result = std::variant<EvaluatedSplitRecords, FilterReason>;

    explicit SplitRecordsEvaluator(
        const std::variant<SplitRecordsEvaluationParameters::BaseParameters,
                           SplitRecordsEvaluationParameters::SplicingParameters> &parameters);

    auto evaluate(SplitRecords &splitRecords, const std::deque<std::string> &referenceIDs) const
        -> Result;

   private:
    std::variant<SplitRecordsEvaluationParameters::BaseParameters,
                 SplitRecordsEvaluationParameters::SplicingParameters>
        parameters;

    static auto evaluateBase(SplitRecords &splitRecords,
                             const SplitRecordsEvaluationParameters::BaseParameters &parameters)
        -> Result;

    static auto evaluateSplicing(
        SplitRecords &splitRecords, const std::deque<std::string> &referenceIDs,
        const SplitRecordsEvaluationParameters::SplicingParameters &parameters) -> Result;

    static void addTagsToRecords(
        SplitRecords &splitRecords,
        const SplitRecordsComplementarityEvaluator::Result &complementarity,
        const SplitRecordsHybridizationEvaluator::Result &hybridization);
};

auto operator<<(std::ostream &ostream, const SplitRecordsEvaluator::FilterReason &reason)
    -> std::ostream &;

}  // namespace pipelines::detect
