#pragma once
#include <map>

// Standard
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <ranges>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

// seqan3
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

// Internal
#include "FigurePlotter.hpp"
#include "HitGroupFailureReason.hpp"
#include "ReadGroup.hpp"
#include "ReadGroupPreprocessor.hpp"
#include "ReadGroupPreprocessorMetrics.hpp"
#include "SegemehlHitGroupConstructor.hpp"

namespace pipelines::detect {

namespace fs = std::filesystem;

struct SegemehlReadGroupPreprocessorConfig {
    SegemehlHitGroupConstructionParameters constructionParams;
    size_t maxPrimaryAlignmentCount;
};

struct SegemehlReadGroupPreprocessor {
    SegemehlReadGroupPreprocessor(SegemehlReadGroupPreprocessorConfig config)
        : collationConfig(config) {}

    [[nodiscard]] constexpr auto getMetrics() const -> const ReadGroupPreprocessorMetrics& {
        return metrics;
    }

    auto operator()(ReadGroup&& readGroup) -> std::vector<ConstructedEvaluationContextVariant> {
        std::unordered_map<int, SegemehlHitGroupConstructor> hitGroupConstructorByTag;
        hitGroupConstructorByTag.reserve(readGroup.size());

        for (auto&& record : std::move(readGroup)) {
            if (record.tags().contains("ch"_tag) || record.tags().contains("SA"_tag)) {
                throw std::runtime_error(
                    "Supplementary alignment tags require supported STAR provenance: " +
                    record.id());
            }
            if (!record.tags().contains("HI"_tag)) {
                throw std::runtime_error("Segemehl alignment lacks required HI tag: " +
                                         record.id());
            }
            const int hitGroupTag = record.tags().get<"HI"_tag>();

            auto [iterator, inserted] = hitGroupConstructorByTag.try_emplace(
                hitGroupTag,
                SegemehlHitGroupConstructor{std::move(record), collationConfig.constructionParams});

            // If it already existed, just insert the moved record
            if (!inserted) {
                iterator->second.insert(std::move(record));
            }
        }

        std::vector<ConstructedEvaluationContextVariant> constructedContexts;
        constructedContexts.reserve(hitGroupConstructorByTag.size());
        for (auto& [tag, constructor] : hitGroupConstructorByTag) {
            constructedContexts.push_back(constructor.getConstructedEvalContext());

            metrics.addContext(constructedContexts.back());
        }

        if (primaryAlignmentCount(constructedContexts) <=
            collationConfig.maxPrimaryAlignmentCount) {
            metrics.finalizeReadGroup();
            return constructedContexts;
        }

        // In case there are more alignments then the cutoff allows for
        std::vector<ConstructedEvaluationContextVariant> errorContextVariants;
        errorContextVariants.reserve(constructedContexts.size());

        for (auto&& contextVariant : constructedContexts) {
            std::visit(
                [&](auto&& context) {
                    using ctx_t = decltype(context);

                    if constexpr (!context.isFailed()) {
                        errorContextVariants.push_back(
                            std::forward<ctx_t>(context).with(HitGroupFailureReason::MULTIMAPPING));
                    } else {
                        errorContextVariants.push_back(std::forward<ctx_t>(context));
                    }
                },
                std::move(contextVariant));
        }

        return errorContextVariants;
    }

   private:
    SegemehlReadGroupPreprocessorConfig collationConfig;
    ReadGroupPreprocessorMetrics metrics{};

    [[nodiscard]] static auto primaryAlignmentCount(
        const std::vector<ConstructedEvaluationContextVariant>& contextVariants) -> size_t {
        std::unordered_map<size_t, size_t> freq;
        freq.reserve(contextVariants.size());
        size_t maxCount = 0;

        for (auto const& contextVariant : contextVariants) {
            size_t mapQ = std::visit(
                [](auto const& ctx) {
                    return ctx.group->getRecordContainer().front().mapping_quality();
                },
                contextVariant);

            size_t& cnt = freq[mapQ];
            ++cnt;
            maxCount = std::max(cnt, maxCount);
        }

        assert(maxCount > 0);
        return maxCount;
    }
};

static_assert(ReadGroupPreprocessor<SegemehlReadGroupPreprocessor>);

}  // namespace pipelines::detect
