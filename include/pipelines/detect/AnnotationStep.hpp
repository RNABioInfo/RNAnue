#pragma once

// Standard
#include <concepts>
#include <memory>
#include <optional>
#include <utility>

// seqan3
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

// Internal
#include "EvaluationContext.hpp"
#include "FeatureAnnotator.hpp"
#include "GenomicFeature.hpp"
#include "GenomicOrientation.hpp"
#include "HitGroup.hpp"
#include "ReadGroupEvaluationParameters.hpp"
#include "SamRecord.hpp"
#include "SplitRecords.hpp"
#include "VariantMeta.hpp"

namespace pipelines::detect {

using namespace dataTypes;

struct AnnotationStepConfig {
    AnnotationStepConfig(std::shared_ptr<const annotation::FeatureAnnotator> featureAnnotator,
                         GenomicOrientation featureOrientation)
        : featureAnnotator(std::move(featureAnnotator)), featureOrientation(featureOrientation) {};

    std::shared_ptr<const annotation::FeatureAnnotator> featureAnnotator;

    GenomicOrientation featureOrientation;

    template <ReadGroupEvaluationParameters::Type EvalParamT>
    [[nodiscard]] static auto makeConfig(const EvalParamT& params) -> AnnotationStepConfig {
        if constexpr (std::same_as<EvalParamT, ReadGroupEvaluationParameters::Base>) {
            return AnnotationStepConfig{params.featureAnnotator, params.featureOrientation};
        } else {
            return AnnotationStepConfig{params.baseParameters.featureAnnotator,
                                        params.baseParameters.featureOrientation};
        }
    }
};

struct AnnotationStepResult {
    std::optional<GenomicFeature> feature;

    void addTags(SamRecord& record) const noexcept {
        if (feature) {
            record.tags()["XF"_tag] = feature->getAnnotationID();
        }
    }
};

struct AnnotationStep {
    using step_result_t = AnnotationStepResult;

    AnnotationStepConfig config;

    template <typename... OtherResults>
    auto operator()(EvaluationContext<HitGroup<SingletonRecord>, OtherResults...>&& ctx) const
        -> EvaluationContext<SingletonHitGroup, OtherResults..., AnnotationStepResult> {
        using ctx_t = decltype(ctx);

        auto recordFeature = config.featureAnnotator->getBestOverlappingFeature(
            ctx.group->getRecordContainer().getRecord(), config.featureOrientation);

        return std::forward_like<ctx_t>(ctx).with(AnnotationStepResult{std::move(recordFeature)});
    }
};

}  // namespace pipelines::detect
