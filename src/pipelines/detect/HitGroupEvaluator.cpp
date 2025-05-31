// #include "HitGroupEvaluator.hpp"

// // Standard
// #include <cassert>
// #include <cstddef>
// #include <cstdint>
// #include <cstdlib>
// #include <type_traits>

// // seqan3
// #include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
// #include <utility>

// // Internal
// #include "CustomSamTags.hpp"
// #include "EvaluatedHitGroup.hpp"
// #include "FailedHitGroup.hpp"
// #include "HitGroup.hpp"
// #include "HitGroupFailureReason.hpp"
// #include "HitGroupFilterReason.hpp"
// #include "LogLevel.hpp"
// #include "Logger.hpp"
// #include "SamRecord.hpp"
// #include "SplitRecords.hpp"
// #include "SplitRecordsComplementarityEvaluator.hpp"
// #include "SplitRecordsEvaluationParameters.hpp"
// #include "SplitRecordsHybridizationEvaluator.hpp"
// #include "SplitRecordsSplicingEvaluator.hpp"

// namespace pipelines::detect::HitGroupEvaluator {

// template <SplitRecordsContainer T>
// auto evaluate(HitGroup<T> &&hitGroup, const HitGroupEvaluation::Parameters::Base &parameters)
//     -> Generator<HitGroupEvaluationResultVariant> {
//     // TODO: Change if including multimeric records
//     if constexpr (not std::is_same_v<T, ChimericRecords>) {
//         co_return EvaluatedHitGroup<T>(hitGroup);
//     }

//     const auto complementarityResult =
//         SplitRecordsComplementarityEvaluator::evaluate(hitGroup.getRecords(), parameters);

//     if (!complementarityResult.has_value()) {
//         co_return FailedHitGroup{std::move(hitGroup),
//                                  HitGroupFailureReason::FAILED_COMPLEMENTARITY};
//     }

//     const auto hybridizationResult =
//         SplitRecordsHybridizationEvaluator::evaluate(hitGroup.getRecords(), parameters);

//     if (!hybridizationResult.has_value()) {
//         co_return FailedHitGroup{std::move(hitGroup),
//         HitGroupFailureReason::FAILED_HYBRIDIZATION};
//     }

//     addTagsToRecords(hitGroup.getRecords(), complementarityResult.value(),
//                      hybridizationResult.value());

//     co_yield HitGroupEvaluation::EvaluatedHitGroup{hitGroup, complementarityResult.value(),
//                                                    hybridizationResult.value()};
// };

// template <SplitRecordsContainer T>
// auto evaluate(HitGroup<T> &&hitGroup, const HitGroupEvaluation::Parameters::Splicing &parameters)
//     -> Generator<HitGroupEvaluationResultVariant> {
//     // TODO: Change if including multimeric records
//     if constexpr (not std::is_same_v<T, ChimericRecords>) {
//         co_return EvaluatedHitGroup<T>(hitGroup);
//     }

//     const auto isSplicing =
//         SplitRecordsSplicingEvaluator::isSplicedSplitRecord(hitGroup.getRecords(), parameters);

//     // TODO: SPLICING Should turn SPLIT Records to a SINGLETON RECORD
//     if (isSplicing) {
//         for (const auto deconstructedHitGroup : hitGroup.getDeconstructedHitGroups()) {
//             co_yield EvaluatedHitGroup<SingletonRecord>(hitGroup);
//         }

//         co_return;
//     }

//     auto baseEvaluate = evaluate(std::move(hitGroup), parameters.baseParameters);
//     for (auto &&result : baseEvaluate) {
//         co_yield result;
//     }
// };

// template <SplitRecordsContainer T>
// void addTagsToRecords(T &splitRecords,
//                       const SplitRecordsComplementarityEvaluator::Result &complementarity,
//                       const SplitRecordsHybridizationEvaluator::Result &hybridization) {
//     assert(splitRecords.size() == 2);  // Currently only two split records are supported

//     size_t index = 0;
//     for (auto &record : splitRecords) {
//         // Complementarity tags
//         const int length = static_cast<int>(complementarity.endPositions.first) -
//                            static_cast<int>(complementarity.beginPositions.first);
//         record.tags()["XL"_tag] = length;
//         record.tags()["XC"_tag] = static_cast<float>(complementarity.complementarity);
//         record.tags()["XR"_tag] = static_cast<float>(complementarity.fraction);
//         record.tags()["XS"_tag] = complementarity.score;

//         // Hybridization tags
//         record.tags()["XE"_tag] = static_cast<float>(hybridization.energy);
//         if (hybridization.crosslinkingResult) {
//             record.tags()["XD"_tag] = hybridization.crosslinkingResult->getDotbracket();
//             record.tags()["XO"_tag] =
//                 static_cast<int32_t>(hybridization.crosslinkingResult->getInterCrosslinkingCount());

//             record.tags()["XA"_tag] =
//                 hybridization.crosslinkingResult->getIntraSequenceCrosslinking(index);

//             record.tags()["XI"_tag] =
//                 hybridization.crosslinkingResult->getInterSequenceCrosslinking(index);
//         }

//         ++index;
//     }
// }

// }  // namespace pipelines::detect::HitGroupEvaluator
