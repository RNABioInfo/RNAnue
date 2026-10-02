#pragma once

#include <deque>
#include <map>
#include <string>
#include <vector>

#include "ReadGroup.hpp"
#include "ReadGroupPreprocessor.hpp"
#include "ReadGroupPreprocessorMetrics.hpp"
#include "SegemehlReadGroupPreprocessor.hpp"
#include "StarAlignmentDecoder.hpp"

namespace pipelines::detect {
class StarReadGroupPreprocessor {
   public:
    StarReadGroupPreprocessor(SegemehlReadGroupPreprocessorConfig config,
                              std::deque<std::string> references, std::vector<size_t> lengths = {})
        : config(config), decoder(std::move(references), std::move(lengths)) {}
    [[nodiscard]] auto getMetrics() const -> const ReadGroupPreprocessorMetrics& { return metrics; }
    auto operator()(ReadGroup&& records) -> std::vector<ConstructedEvaluationContextVariant> {
        auto hypotheses = decoder.decode(records, config.constructionParams.excludeSoftClipping);
        std::vector<ConstructedEvaluationContextVariant> contexts;

        for (auto& hypothesis : hypotheses) {
            if (hypothesis.fragments.empty()) {
                continue;
            }
            std::optional<HitGroupFailureReason> failure;

            if (hypothesis.paired) {
                failure = HitGroupFailureReason::UNSUPPORTED_PAIRED;
            } else if (hypothesis.fragments.size() > 2) {
                failure = HitGroupFailureReason::UNSUPPORTED_MULTISEGMENT;
            } else {
                for (const auto& fragment : hypothesis.fragments) {
                    if (fragment.mapping_quality() < config.constructionParams.minimumMapQuality) {
                        failure = HitGroupFailureReason::MAPPING_QUALITY;
                    } else if (fragment.sequence().size() <
                               config.constructionParams.minimumFragmentLength) {
                        failure = HitGroupFailureReason::FRAGMENT_LENGTH;
                    }
                }
            }

            const auto append = [&](auto group) {
                if (failure) {
                    contexts.emplace_back(
                        EvaluationContext{std::move(group), std::make_tuple(*failure)});
                } else {
                    contexts.emplace_back(EvaluationContext{std::move(group)});
                }
            };

            auto& fragments = hypothesis.fragments;

            if (fragments.size() == 1) {
                append(std::make_unique<SingletonHitGroup>(
                    SingletonRecord{std::move(fragments.front())}, hypothesis.hitIndex));
            } else if (fragments.size() == 2) {
                append(std::make_unique<ChimericHitGroup>(
                    ChimericRecords{std::move(fragments.front()), std::move(fragments.back())},
                    hypothesis.hitIndex));
            } else {
                append(std::make_unique<MultimericHitGroup>(MultimericRecords{std::move(fragments)},
                                                            hypothesis.hitIndex));
            }
        }

        // Match the existing MAPQ-group rule, counting hypotheses rather than BAM rows.
        std::map<size_t, size_t> counts;
        size_t maximum = 0;

        for (const auto& context : contexts) {
            const auto quality = std::visit(
                [](const auto& ctx) {
                    return ctx.group->getRecordContainer().front().mapping_quality();
                },
                context);
            maximum = std::max(maximum, ++counts[quality]);
        }

        if (maximum > config.maxPrimaryAlignmentCount) {
            std::vector<ConstructedEvaluationContextVariant> failed;
            failed.reserve(contexts.size());

            for (auto& context : contexts) {
                failed.emplace_back(std::visit(
                    [](auto&& ctx) -> ConstructedEvaluationContextVariant {
                        if constexpr (std::remove_cvref_t<decltype(ctx)>::isFailed()) {
                            return std::forward<decltype(ctx)>(ctx);
                        } else {
                            return std::forward<decltype(ctx)>(ctx).with(
                                HitGroupFailureReason::MULTIMAPPING);
                        }
                    },
                    std::move(context)));
            }
            contexts.swap(failed);
        }

        for (auto& context : contexts) {
            metrics.addContext(context);
        }

        metrics.finalizeReadGroup();
        return contexts;
    }

   private:
    SegemehlReadGroupPreprocessorConfig config;
    align::StarAlignmentDecoder decoder;
    ReadGroupPreprocessorMetrics metrics;
};
static_assert(ReadGroupPreprocessor<StarReadGroupPreprocessor>);
}  // namespace pipelines::detect
