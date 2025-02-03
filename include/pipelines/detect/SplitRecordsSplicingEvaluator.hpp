#pragma once

// Standard
#include <cstddef>
#include <utility>
#include <vector>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "SortedGenomicRegionPair.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"

using namespace dataTypes;
using namespace annotation;

class SplitRecordsSplicingEvaluator {
   public:
    SplitRecordsSplicingEvaluator() = delete;
    ~SplitRecordsSplicingEvaluator() = delete;
    SplitRecordsSplicingEvaluator(const SplitRecordsSplicingEvaluator &) = delete;
    auto operator=(const SplitRecordsSplicingEvaluator &)
        -> SplitRecordsSplicingEvaluator & = delete;
    SplitRecordsSplicingEvaluator(SplitRecordsSplicingEvaluator &&) = delete;
    auto operator=(SplitRecordsSplicingEvaluator &&) -> SplitRecordsSplicingEvaluator & = delete;

    static auto isSplicedSplitRecord(
        const SplitRecords &splitRecords,
        const SplitRecordsEvaluationParameters::SplicingParameters &parameters) -> bool;

   private:
    using BoundingRegion = GenomicRegion;
    using BoundingRegions = std::pair<BoundingRegion, BoundingRegion>;
    using FeaturePair = std::pair<GenomicFeature, GenomicFeature>;
    using FeaturePairs = std::vector<FeaturePair>;

    static auto getGroupedFeaturesAtBoundingRegions(
        const BoundingRegions &boundingRegions,
        const SplitRecordsEvaluationParameters::SplicingParameters &parameters) -> FeaturePairs;

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
    static auto groupedFeaturePairsEncloseAdditonalFeatureFromGroup(
        const FeaturePairs &featurePairs,
        const SplitRecordsEvaluationParameters::SplicingParameters &parameters) -> bool;
};
