#pragma once

// Standard
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <utility>
#include <variant>
#include <vector>

// Internal
#include "EvaluationContext.hpp"
#include "FeatureAnnotator.hpp"
#include "GenomicFeature.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "HitGroup.hpp"
#include "OneOf.hpp"
#include "ReadGroupEvaluationParameters.hpp"
#include "SamRecord.hpp"
#include "SortedGenomicRegionPair.hpp"
#include "SplitRecords.hpp"

namespace pipelines::detect {

using namespace dataTypes;

struct SplicingEvaluationStepConfig {
    SplicingEvaluationStepConfig(
        std::shared_ptr<const annotation::FeatureAnnotator> featureAnnotator,
        GenomicOrientation featureOrientation, int32_t splicingTolerance,
        bool removeAlternativeSplicing)
        : featureAnnotator(std::move(featureAnnotator)),
          featureOrientation(featureOrientation),
          splicingTolerance(splicingTolerance),
          removeAlternativeSplicing(removeAlternativeSplicing) {};

    std::shared_ptr<const annotation::FeatureAnnotator> featureAnnotator;
    GenomicOrientation featureOrientation;
    int32_t splicingTolerance;
    bool removeAlternativeSplicing;

    [[nodiscard]] static auto makeConfig(const ReadGroupEvaluationParameters::Splicing &params)
        -> SplicingEvaluationStepConfig {
        return SplicingEvaluationStepConfig{
            params.baseParameters.featureAnnotator, params.baseParameters.featureOrientation,
            params.splicingTolerance, params.allowAlternativeSplicing};
    };
};

using FeaturePair = std::pair<GenomicFeature, GenomicFeature>;
using FeaturePairs = std::vector<FeaturePair>;

struct SplicingEvaluationStepResult {
    std::optional<FeaturePairs> splicingFeatures;

    [[nodiscard]] constexpr auto isSpliced() const noexcept -> bool {
        return splicingFeatures.has_value();
    }

    void addTags(SamRecord & /*unused */) const {}
};

class SplicingEvaluationStep {
   public:
    using step_result_t = SplicingEvaluationStepResult;

    template <typename... OtherResults>
    using step_return_t = std::variant<
        std::vector<EvaluationContext<SingletonHitGroup, OtherResults..., step_result_t>>,
        EvaluationContext<ChimericHitGroup, OtherResults..., step_result_t>>;

    SplicingEvaluationStepConfig config;

    template <typename... OtherResults>
    auto operator()(EvaluationContext<ChimericHitGroup, OtherResults...> &&ctx) const
        -> step_return_t<OtherResults...> {
        using ctx_t = decltype(ctx);

        const std::optional<GenomicRegion> regionOne =
            GenomicRegion::fromSamRecord(ctx.group->getRecordContainer().first());
        const std::optional<GenomicRegion> regionTwo =
            GenomicRegion::fromSamRecord(ctx.group->getRecordContainer().second());

        if (!regionOne.has_value() || !regionTwo.has_value()) {
            return {std::forward_like<ctx_t>(ctx).with(SplicingEvaluationStepResult{})};
        }

        const auto boundingRegions =
            getBoundingRegionsForRecords({*regionOne, *regionTwo}, config.splicingTolerance);

        auto featurePairs = getGroupedFeaturesAtBoundingRegions(boundingRegions);

        if (featurePairs.empty()) {
            return {std::forward_like<ctx_t>(ctx).with(SplicingEvaluationStepResult{})};
        }

        // There are feature pairs and it doesn't matter if there are exons in between
        if (config.removeAlternativeSplicing ||
            !groupedFeaturePairsEncloseAdditonalFeatureFromGroup(featurePairs)) {
            return std::forward_like<ctx_t>(ctx)
                .with(SplicingEvaluationStepResult{featurePairs})
                .getDeconstructedContexts();
        }

        return std::forward_like<ctx_t>(ctx).with(SplicingEvaluationStepResult{});
    };

    template <typename T>
        requires one_of<T, ChimericRecords, MultimericRecords>
    static auto isSplicedSplitRecord(const T &splitRecords,
                                     const ReadGroupEvaluationParameters::Splicing &parameters)
        -> bool;

   private:
    using BoundingRegion = GenomicRegion;
    using BoundingRegions = std::pair<BoundingRegion, BoundingRegion>;

    [[nodiscard]] auto getGroupedFeaturesAtBoundingRegions(
        const BoundingRegions &boundingRegions) const -> FeaturePairs;

    /** @brief Get the bounding regions for two records.
     *
     * This function calculates the bounding regions for two records based on the reference
     * positions of the records. No tolerance means that only the exact end and start positions are
     * returned.
     *
     * @param record1 The first record.
     * @param record2 The second record.
     * @param tolerance The tolerance for the bounding regions.
     * @return The bounding regions for the two records.
     */
    static auto getBoundingRegionsForRecords(const SortedGenomicRegionPair &regionPair,
                                             const size_t &tolerance) -> BoundingRegions;

    /**
     * @brief Checks if grouped feature pairs enclose an additional feature from the same group.
     *
     * This function iterates through each feature pair and constructs a region bounded by
     * the end of the first feature and the start of the second. It then checks if there is
     * another feature in this region belonging to the same group as the first feature.
     *
     * @param featurePairs The collection of feature pairs to examine.
     * @param parameters The splicing parameters that provide the feature annotator and orientation.
     * @return True if any bounding region encloses a feature of the same group as the first feature
     * in the pair.
     */
    [[nodiscard]] auto groupedFeaturePairsEncloseAdditonalFeatureFromGroup(
        const FeaturePairs &featurePairs) const -> bool;
};

}  // namespace pipelines::detect
