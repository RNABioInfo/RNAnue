#include "SplitRecordsEvaluator.hpp"

// Standard
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <ostream>
#include <variant>

// seqan3
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

// Internal
#include "CustomSamTags.hpp"
#include "Logger.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsComplementarityEvaluator.hpp"
#include "SplitRecordsEvaluationParameters.hpp"
#include "SplitRecordsHybridizationEvaluator.hpp"
#include "SplitRecordsSplicingEvaluator.hpp"

namespace pipelines::detect {

SplitRecordsEvaluator::SplitRecordsEvaluator(
    const std::variant<SplitRecordsEvaluationParameters::BaseParameters,
                       SplitRecordsEvaluationParameters::SplicingParameters> &parameters)
    : parameters(parameters) {}

auto SplitRecordsEvaluator::evaluate(SplitRecords &splitRecords) const
    -> SplitRecordsEvaluator::Result {
    if (splitRecords.size() != 2) {
        Logger::log<LogLevel::DEBUG>("Currently only two split records are supported!");
        return SplitRecordsEvaluator::FilterReason::NO_SPLIT_READ;
    }

    if (!splitRecords[0].reference_position().has_value() ||
        !splitRecords[1].reference_position().has_value()) [[unlikely]] {
        Logger::log<LogLevel::WARNING>("Could not determine reference position of split records: ",
                                       splitRecords[0].id(), " or ", splitRecords[1].id());
        return SplitRecordsEvaluator::FilterReason::UNMAPPED;
    }

    if (std::holds_alternative<SplitRecordsEvaluationParameters::BaseParameters>(parameters)) {
        return evaluateBase(splitRecords,
                            std::get<SplitRecordsEvaluationParameters::BaseParameters>(parameters));
    }

    return evaluateSplicing(
        splitRecords, std::get<SplitRecordsEvaluationParameters::SplicingParameters>(parameters));
}

auto SplitRecordsEvaluator::evaluateBase(
    SplitRecords &splitRecords, const SplitRecordsEvaluationParameters::BaseParameters &parameters)
    -> SplitRecordsEvaluator::Result {
    const auto complementarityResult =
        SplitRecordsComplementarityEvaluator::evaluate(splitRecords, parameters);

    if (!complementarityResult.has_value()) {
        return SplitRecordsEvaluator::FilterReason::COMPLEMENTARITY;
    }

    const auto hybridizationResult =
        SplitRecordsHybridizationEvaluator::evaluate(splitRecords, parameters);

    if (!hybridizationResult.has_value()) {
        return SplitRecordsEvaluator::FilterReason::HYBRIDIZATION;
    }

    addTagsToRecords(splitRecords, complementarityResult.value(), hybridizationResult.value());

    return SplitRecordsEvaluator::EvaluatedSplitRecords{
        .splitRecords = splitRecords,
        .complementarityResult = complementarityResult.value(),
        .hybridizationResult = hybridizationResult.value()};
}

auto SplitRecordsEvaluator::evaluateSplicing(
    SplitRecords &splitRecords,
    const SplitRecordsEvaluationParameters::SplicingParameters &parameters)
    -> SplitRecordsEvaluator::Result {
    const auto isSplicing =
        SplitRecordsSplicingEvaluator::isSplicedSplitRecord(splitRecords, parameters);

    if (isSplicing) {
        return SplitRecordsEvaluator::FilterReason::SPLICING;
    }

    return evaluateBase(splitRecords, parameters.baseParameters);
}

void SplitRecordsEvaluator::addTagsToRecords(
    SplitRecords &splitRecords, const SplitRecordsComplementarityEvaluator::Result &complementarity,
    const SplitRecordsHybridizationEvaluator::Result &hybridization) {
    assert(splitRecords.size() == 2);  // Currently only two split records are supported

    size_t index = 0;
    for (auto &record : splitRecords) {
        // Complementarity tags
        const int length = static_cast<int>(complementarity.endPositions.first) -
                           static_cast<int>(complementarity.beginPositions.first);
        record.tags()["XL"_tag] = length;
        record.tags()["XC"_tag] = static_cast<float>(complementarity.complementarity);
        record.tags()["XR"_tag] = static_cast<float>(complementarity.fraction);
        record.tags()["XS"_tag] = complementarity.score;

        // Hybridization tags
        record.tags()["XE"_tag] = static_cast<float>(hybridization.energy);
        if (hybridization.crosslinkingResult) {
            record.tags()["XD"_tag] = hybridization.crosslinkingResult->getDotbracket();
            record.tags()["XO"_tag] =
                static_cast<int32_t>(hybridization.crosslinkingResult->getTotalCrosslinkingCount());

            record.tags()["XA"_tag] =
                hybridization.crosslinkingResult->getIntraSequenceCrosslinking(index);

            record.tags()["XI"_tag] =
                hybridization.crosslinkingResult->getInterSequenceCrosslinking(index);
        }

        ++index;
    }
}

// TODO Implement overall score calculation instead of comparing individual scores
auto SplitRecordsEvaluator::EvaluatedSplitRecords::operator<(
    const EvaluatedSplitRecords &other) const -> bool {
    constexpr double DIFF_THRESHOLD = 0.05;

    // Calculate if fraction difference is more or equal than threshold
    if (std::abs(complementarityResult.fraction - other.complementarityResult.fraction) >=
        DIFF_THRESHOLD) {
        return complementarityResult.fraction < other.complementarityResult.fraction;
    }

    // Calculate if complementarity difference is more or equal than threshold
    if (std::abs(complementarityResult.complementarity -
                 other.complementarityResult.complementarity) >= DIFF_THRESHOLD) {
        return complementarityResult.complementarity < other.complementarityResult.complementarity;
    }

    // Calculate if energy difference is more or equal than threshold
    if (std::abs(hybridizationResult.energy - other.hybridizationResult.energy) >= DIFF_THRESHOLD) {
        return hybridizationResult.energy < other.hybridizationResult.energy;
    }

    // If both have crosslinking results prefer the one with higher crosslinking score
    if (hybridizationResult.crosslinkingResult.has_value() &&
        other.hybridizationResult.crosslinkingResult.has_value()) {
        return hybridizationResult.crosslinkingResult->getTotalCrosslinkingCount() <
               other.hybridizationResult.crosslinkingResult->getTotalCrosslinkingCount();
    }

    return false;
}

auto SplitRecordsEvaluator::EvaluatedSplitRecords::operator>(
    const EvaluatedSplitRecords &other) const -> bool {
    return !(*this < other);
}

auto operator<<(std::ostream &ostream, const SplitRecordsEvaluator::FilterReason &reason)
    -> std::ostream & {
    switch (reason) {
        case SplitRecordsEvaluator::FilterReason::NO_SPLIT_READ:
            ostream << "No split read";
            break;
        case SplitRecordsEvaluator::FilterReason::UNMAPPED:
            ostream << "Unmapped";
            break;
        case SplitRecordsEvaluator::FilterReason::SPLICING:
            ostream << "Splicing";
            break;
        case SplitRecordsEvaluator::FilterReason::COMPLEMENTARITY:
            ostream << "Complementarity";
            break;
        case SplitRecordsEvaluator::FilterReason::HYBRIDIZATION:
            ostream << "Hybridization";
            break;
    }

    return ostream;
}

}  // namespace pipelines::detect
