#pragma once

// Standard
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <format>
#include <iterator>
#include <limits>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Seqan3
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/utility/views/slice.hpp>

// Internal
#include "AlignParameters.hpp"
#include "GenomicFeature.hpp"
#include "GenomicFeatureGroup.hpp"
#include "GenomicRegion.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "MaskedFeatureCluster.hpp"
#include "ReferenceGenome.hpp"
#include "seqan3/alignment/configuration/align_config_min_score.hpp"
#include "seqan3/alignment/configuration/align_config_output.hpp"

namespace pipelines::align {

class GeneCopyMasker {
   public:
    using FeatureID = std::string;

    struct Parameters {
        explicit Parameters(const AlignParameters& params)
            : maskMultiCopyGenes(params.maskMultiCopyGenes),
              minMultiCopyIdentity(params.minMultiCopyIdentity) {}

        bool maskMultiCopyGenes{};
        double minMultiCopyIdentity{};
    };

    struct Result {
        ReferenceGenome maskedGenome;
        std::vector<MaskedFeatureCluster> maskedClusters;
    };

    explicit GeneCopyMasker(const Parameters& params, ParentIDToFeatureGroupMap&& featureGroups,
                            ReferenceGenome&& referenceGenomeParam)
        : parameters{params},
          featureGroups{std::move(featureGroups)},
          referenceGenome{std::move(referenceGenomeParam)},
          maskedGenome(this->referenceGenome) {}

    [[nodiscard]] auto process() -> Result {
        if (!parameters.maskMultiCopyGenes) {
            return Result{.maskedGenome = std::move(maskedGenome),
                          .maskedClusters = std::move(clusters)};
        }

        Logger::log(
            std::format("Starting gene copy masking for {} features.", featureGroups.size()));
        const auto alignmentConfig = seqan3::align_cfg::method_global{} |
                                     seqan3::align_cfg::edit_scheme |
                                     seqan3::align_cfg::output_score{};

        clusters.reserve(featureGroups.size());

        size_t numMasked = 0;

        for (auto& group : std::views::values(featureGroups)) {
            if (maskGroupIfCopy(group, alignmentConfig)) {
                ++numMasked;
            };
        }

        Logger::log(std::format("Masked {} features.", numMasked));

        return Result{.maskedGenome = std::move(maskedGenome),
                      .maskedClusters = std::move(clusters)};
    }

   private:
    static constexpr auto kMinIdentityBound = 0.0;
    static constexpr auto kMaxIdentityBound = 1.0;
    static constexpr auto kFullIdentity = 1.0;
    static constexpr auto kZeroEdits = std::size_t{0};
    static constexpr auto kMaskNucleotide = 'N'_dna5;

    Parameters parameters;
    ParentIDToFeatureGroupMap featureGroups;
    const ReferenceGenome referenceGenome;
    ReferenceGenome maskedGenome;

    std::vector<MaskedFeatureCluster> clusters;
    std::unordered_map<std::size_t, std::vector<std::size_t>> clusterIndicesByLength;

    template <typename TAlignmentConfig>
    [[nodiscard]] auto maskGroupIfCopy(dataTypes::GenomicFeatureGroup& group,
                                       TAlignmentConfig const& alignmentConfig) -> bool {
        const auto& rootFeature = group.getRoot().feature;
        const auto rootRegion = rootFeature.getGenomicRegion();

        const auto querySequence = extractRegionSequenceSpan(rootRegion);
        if (querySequence.empty()) {
            Logger::log<LogLevel::WARNING>(
                std::format("Found empty query sequence for {}", rootFeature.getAnnotationID()));
            return false;
        }

        const auto currentLength = querySequence.size();

        const auto [minCandidateLength, maxCandidateLength] =
            allowedLengthRange(currentLength, parameters.minMultiCopyIdentity);

        if (const auto matchedClusterIndex = findMatchingClusterIndex(
                querySequence, minCandidateLength, maxCandidateLength, alignmentConfig);
            matchedClusterIndex.has_value()) {
            Logger::log<LogLevel::DEBUG>(
                std::format("Found match at index: {}, for feature group: {}", *matchedClusterIndex,
                            group.getRoot().feature.getAnnotationID()));
            registerClusterMatch(clusters[*matchedClusterIndex], group);
            return true;
        }

        addNewCluster(seqan3::dna5_vector{querySequence.begin(), querySequence.end()}, group);

        return false;
    };

    template <typename TAlignmentConfig>
    [[nodiscard]] auto findMatchingClusterIndex(std::span<const seqan3::dna5> currentSequence,
                                                std::size_t minCandidateLength,
                                                std::size_t maxCandidateLength,
                                                TAlignmentConfig const& alignmentConfig)
        -> std::optional<std::size_t> {
        for (const auto& [bucketLength, bucketClusterIndices] : clusterIndicesByLength) {
            if (bucketLength < minCandidateLength || bucketLength > maxCandidateLength) {
                continue;
            }

            for (const auto clusterIndex : bucketClusterIndices) {
                if (clusterMatches(clusters[clusterIndex], currentSequence, alignmentConfig)) {
                    return clusterIndex;
                }
            }
        }

        return std::nullopt;
    }

    template <typename TAlignmentConfig>
    [[nodiscard]] auto clusterMatches(MaskedFeatureCluster const& cluster,
                                      std::span<const seqan3::dna5> currentSequence,
                                      TAlignmentConfig const& alignmentConfig) const -> bool {
        const auto currentLength = currentSequence.size();
        const auto candidateLength = cluster.baseSequence.size();

        const auto requiredEditsByLength = (candidateLength > currentLength)
                                               ? (candidateLength - currentLength)
                                               : (currentLength - candidateLength);

        const auto allowedEdits =
            maxAllowedEdits(currentLength, candidateLength, parameters.minMultiCopyIdentity);

        Logger::log<LogLevel::DEBUG>(std::format(
            "Computing edit distance with settings. requiredEditsByLength: {}, allowedEdits: {}",
            requiredEditsByLength, allowedEdits));

        if (requiredEditsByLength > allowedEdits) {
            return false;
        }

        if (std::ranges::equal(cluster.baseSequence, currentSequence)) {
            Logger::log<LogLevel::DEBUG>("Got exact match for sequences.");
            return true;
        }

        const auto editCountOpt = computeEditDistance(cluster.baseSequence, currentSequence,
                                                      alignmentConfig, allowedEdits);

        if (!editCountOpt.has_value()) {
            return false;
        }

        const auto identity = editDistanceToIdentity(*editCountOpt, currentLength, candidateLength);

        Logger::log<LogLevel::DEBUG>(std::format("Gene identity is {}, min identity is {}",
                                                 identity, parameters.minMultiCopyIdentity));

        return identity >= parameters.minMultiCopyIdentity;
    }

    auto registerClusterMatch(MaskedFeatureCluster& cluster,
                              GenomicFeatureGroup const& featureGroup) -> void {
        cluster.subFeatureGroups.emplace_back(featureGroup);

        maskRegion(featureGroup.getRoot().feature.getGenomicRegion());
    };

    auto addNewCluster(seqan3::dna5_vector&& representativeSequence,
                       GenomicFeatureGroup const& featureGroup) -> void;

    [[nodiscard]] auto extractRegionSequenceSpan(const GenomicRegion& region) const
        -> seqan3::dna5_vector;

    auto maskRegion(const GenomicRegion& region) -> void;

    [[nodiscard]] static auto allowedLengthRange(std::size_t length, double minIdentity)
        -> std::pair<std::size_t, std::size_t>;

    [[nodiscard]] static auto maxAllowedEdits(std::size_t lengthFirst, std::size_t lengthSecond,
                                              double minIdentity) -> std::size_t;

    [[nodiscard]] static auto editDistanceToIdentity(std::size_t editCount, std::size_t lengthFirst,
                                                     std::size_t lengthSecond) -> double;

    template <typename TAlignmentConfig, typename TSeqA, typename TSeqB>
    [[nodiscard]] static auto computeEditDistance(TSeqA const& sequenceFirst,
                                                  TSeqB const& sequenceSecond,
                                                  TAlignmentConfig const& alignmentConfig,
                                                  int32_t allowedEdits)
        -> std::optional<std::size_t> {
        static constexpr auto kInfinityScore = std::numeric_limits<int32_t>::max();

        // Clamp allowedEdits to a sane non-negative range.
        if (allowedEdits < 0) {
            return std::nullopt;
        }

        const auto minScoreSigned = [&]() -> int32_t {
            const auto clipped = std::min(allowedEdits, kInfinityScore);
            return -clipped;  // edit distance score is negative edits
        }();

        const auto alignmentConfigCopy =
            alignmentConfig | seqan3::align_cfg::min_score{minScoreSigned};

        auto results =
            seqan3::align_pairwise(std::tie(sequenceFirst, sequenceSecond), alignmentConfigCopy);

        const auto iter = std::ranges::begin(results);
        if (iter == std::ranges::end(results)) {
            return std::nullopt;
        }

        const auto score = iter->score();

        // SeqAn: if bound is violated, score is "infinity" == numeric_limits::max().
        if (score == kInfinityScore) {
            return std::nullopt;
        }

        // With edit_scheme, best is 0, and worse scores are negative.
        if (score > 0) {  // should not happen for edit distance; fail closed
            return std::nullopt;
        }

        return score == 0 ? kZeroEdits : static_cast<std::size_t>(-static_cast<int64_t>(score));
    }
};

}  // namespace pipelines::align
