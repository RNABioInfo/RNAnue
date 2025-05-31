// #pragma once

// // Internal
// #include <utility>

// #include "EvaluatedHitGroup.hpp"
// #include "FailedHitGroup.hpp"
// #include "HitGroup.hpp"
// #include "HitGroupConstructor.hpp"
// #include "HitGroupEvaluationResult.hpp"
// #include "HitGroupFilterReason.hpp"
// #include "SplitRecords.hpp"
// #include "SplitRecordsComplementarityEvaluator.hpp"
// #include "SplitRecordsEvaluationParameters.hpp"
// #include "SplitRecordsHybridizationEvaluator.hpp"
// #include "SplitRecordsSplicingEvaluator.hpp"
// #include "VariantUnion.hpp"

// using namespace seqan3::literals;

// namespace pipelines::detect {

// using namespace dataTypes;

// using HitGroupEvaluationResultVariant =
//     VariantUnion<HitGroupEvaluation::EvaluatedHitGroupVariant, FailedHitGroupVariant>;

// namespace HitGroupEvaluator {
// using FilterReason = HitGroupFilterReason;

// namespace {
// template <SplitRecordsContainer T>
// void addTagsToRecords(T &splitRecords,
//                       const HitGroupEvaluation::ComplementarityResult &complementarity,
//                       const HitGroupEvaluation::HybridizationResult &hybridization);
// }  // namespace

// template <SplitRecordsContainer T>
// [[nodiscard]] auto evaluate(HitGroup<T> &&hitGroup,
//                             const HitGroupEvaluation::Parameters::Base &parameters)
//     -> Generator<HitGroupEvaluationResultVariant>;

// template <SplitRecordsContainer T>
// [[nodiscard]] auto evaluate(HitGroup<T> &&hitGroup,
//                             const HitGroupEvaluation::Parameters::Splicing &parameters)
//     -> Generator<HitGroupEvaluationResultVariant>;

// };  // namespace HitGroupEvaluator

// }  // namespace pipelines::detect
