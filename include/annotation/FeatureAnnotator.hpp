#pragma once

// Boost
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard
#include <cstddef>
#include <filesystem>
#include <iterator>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicFeatureTreeMerger.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "IITree.hpp"
#include "SamRecord.hpp"

namespace annotation {

using namespace dataTypes;

/// Mapping from a reference name to its index.
using ReferenceIDIndex = int;
using ReferenceIDToIndexMap = std::unordered_map<std::string, int>;

/// Map from reference index to the interval tree of features.
using FeatureTreeMap = std::unordered_map<ReferenceIDIndex, IITree<int, GenomicFeature>>;

namespace fs = std::filesystem;

/**
 * @brief Loads (and checks consistency for) the mapping from reference sequence
 *        identifiers to indices from a list of SAM file paths.
 *
 * @param samFilePaths Vector of paths to SAM files.
 * @return ReferenceIDToIndexMap Mapping from reference ID to index.
 */
[[nodiscard]] auto loadReferenceIDToIndexMap(const std::vector<fs::path>& samFilePaths)
    -> ReferenceIDToIndexMap;

/**
 * @brief Class for annotating and managing genomic features.
 */
class FeatureAnnotator {
   public:
    /**
     * @brief Construct a new Feature Annotator object using a feature file.
     *
     * @param featureFilePath Path to the feature file.
     * @param referenceIDToIndex Pre-built mapping from reference IDs to indices.
     * @param includedFeatures Set of feature names to include.
     * @param featureIDFlag Feature identifier flag.
     */
    FeatureAnnotator(fs::path& featureFilePath, const ReferenceIDToIndexMap& referenceIDToIndex,
                     const std::unordered_set<std::string>& includedFeatures,
                     const std::string& featureIDFlag) noexcept;
    /**
     * @brief Construct a new Feature Annotator object using a feature file.
     *
     * @param featureFilePath Path to the feature file.
     * @param referenceIDToIndex Pre-built mapping from reference IDs to indices.
     * @param includedFeatures Set of feature names to include.
     */
    FeatureAnnotator(fs::path& featureFilePath, const ReferenceIDToIndexMap& referenceIDToIndex,
                     const std::unordered_set<std::string>& includedFeatures) noexcept;
    /**
     * @brief Construct a Feature Annotator object using an in-memory feature map.
     *
     * @param featureMap Map from reference index to vector of GenomicFeature.
     */
    explicit FeatureAnnotator(const FeatureMap& featureMap) noexcept;

    explicit FeatureAnnotator() = default;

    FeatureAnnotator(const FeatureAnnotator&) = default;
    FeatureAnnotator(FeatureAnnotator&&) = default;
    auto operator=(const FeatureAnnotator&) -> FeatureAnnotator& = default;
    auto operator=(FeatureAnnotator&&) -> FeatureAnnotator& = delete;
    ~FeatureAnnotator() = default;

    /// Forward iterator helper for a collection of overlapping features.
    class Results;
    struct MergeInsertResult;

    /**
     * @brief Returns the total feature count.
     *
     * @return size_t Total number of features.
     */
    [[nodiscard]] auto featureCount() const noexcept -> size_t;

    /**
     * @brief Inserts an index covering a given genomic region.
     *
     * @param region Genomic region to index.
     * @return std::string A newly generated feature identifier.
     */
    auto insertIndex(const GenomicRegion& region) -> std::string;

    /**
     * @brief Returns a lazy iterable of features overlapping a given region.
     *
     * @param region Genomic region of interest.
     * @param orientation The genomic orientation (SAME, OPPOSITE, or BOTH).
     * @return Results Iterable results over overlapping features.
     */
    [[nodiscard]] auto overlappingFeatureIt(const GenomicRegion& region,
                                            GenomicOrientation orientation) const -> Results;

    /**
     * @brief Returns vector of features overlapping a given region.
     *
     * @param region Genomic region.
     * @param orientation The genomic orientation.
     * @return std::vector<GenomicFeature> Overlapping features.
     */
    [[nodiscard]] auto getOverlappingFeatures(const GenomicRegion& region,
                                              GenomicOrientation orientation) const
        -> std::vector<GenomicFeature>;
    /**
     * @brief Returns the best overlapping feature for a region.
     *
     * @param region Genomic region.
     * @param orientation The genomic orientation.
     * @return std::optional<GenomicFeature> The best (highest overlap) feature, if any.
     */
    [[nodiscard]] auto getBestOverlappingFeature(const GenomicRegion& region,
                                                 GenomicOrientation orientation) const
        -> std::optional<GenomicFeature>;
    /**
     * @brief Returns the best overlapping feature for a SAM record.
     *
     * @param record SAM record.
     * @param orientation The genomic orientation.
     * @return std::optional<GenomicFeature> The best (highest overlap) feature, if any.
     */
    [[nodiscard]] auto getBestOverlappingFeature(const SamRecord& record,
                                                 GenomicOrientation orientation) const
        -> std::optional<GenomicFeature>;

    /**
     * @brief Returns the best overlapping feature for a SAM record.
     *
     * @param record SAM record.
     * @param orientation The preferred genomic orientation.
     * @return std::optional<GenomicFeature> The best (highest overlap) feature, if any.
     */
    [[nodiscard]] auto getBestOverlappingFeatureWithPreferredOrientation(
        const GenomicRegion& region, GenomicOrientation orientation) const
        -> std::optional<GenomicFeature>;

    /**
     * @brief Getter for the internal feature tree map.
     *
     * @return const FeatureTreeMap& Constant reference to the feature tree map.
     */
    [[nodiscard]] auto getFeatureTreeMap() const -> const FeatureTreeMap&;

    /**
     * @brief Merges all overlapping features in the annotation trees.
     *
     * @param mergingParameters Parameters controlling the merging behavior.
     */
    void mergeAllOverlappingFeatures(FeatureMergingParameters mergingParameters);

    /**
     * @brief Outputs all features to std::cout for debugging.
     */
    void debugOutputAllFeatures() const;

   private:
    FeatureTreeMap featureTreeMap;

    static auto buildFeatureTreeMap(const fs::path& featureFilePath,
                                    const ReferenceIDToIndexMap& referenceIDToIndex,
                                    const std::vector<std::string>& includedFeatures,
                                    const std::optional<std::string>& featureIDFlag)
        -> FeatureTreeMap;

    /**
     * @brief Parses and builds the feature tree map from a feature file.
     *
     * @param featureFilePath Path to the feature file.
     * @param referenceIDToIndex Mapping from reference ID to index.
     * @param includedFeatures Set of feature names to include.
     * @param featureIDFlag Optional flag to identify features.
     * @return FeatureTreeMap The constructed tree map.
     */
    static auto buildFeatureTreeMap(const fs::path& featureFilePath,
                                    const ReferenceIDToIndexMap& referenceIDToIndex,
                                    const std::unordered_set<std::string>& includedFeatures,
                                    const std::optional<std::string>& featureIDFlag)
        -> FeatureTreeMap;

    /**
     * @brief Builds a feature tree map from an in-memory feature map.
     *
     * @param featureMap Map from reference index to vector of GenomicFeature.
     * @return FeatureTreeMap The constructed tree map.
     */
    static auto buildFeatureTreeMap(const FeatureMap& featureMap) -> FeatureTreeMap;
};

/// Forward iterator over a collection of overlapping genomic features.
class FeatureAnnotator::Results {
   public:
    /**
     * @brief Construct a new Results object.
     *
     * @param tree Pointer to the IITree containing the features (may be nullptr if none).
     * @param indices Vector of indices corresponding to overlapping features.
     * @param strand Optional filter strand.
     */
    Results(const IITree<int, GenomicFeature>* tree, const std::vector<size_t>& indices,
            std::optional<GenomicStrand> strand);

    Results() = delete;

    struct Iterator;

    /**
     * @brief Returns iterator to the beginning.
     *
     * @return Iterator Starting iterator.
     */
    [[nodiscard]] auto begin() const -> Iterator;
    /**
     * @brief Returns iterator to the end.
     *
     * @return Iterator Ending iterator.
     */
    [[nodiscard]] auto end() const -> Iterator;

   private:
    const IITree<int, GenomicFeature>* tree;
    std::vector<size_t> indices;
    std::optional<GenomicStrand> strand;
};

/// Iterator class for Results.
struct FeatureAnnotator::Results::Iterator {
    using value_type = GenomicFeature;
    using difference_type = std::ptrdiff_t;
    using reference = const GenomicFeature&;
    using pointer = const GenomicFeature*;
    using iterator_category = std::forward_iterator_tag;

    /**
     * @brief Construct an Iterator.
     *
     * @param tree Pointer to the underlying IITree.
     * @param indices Vector of indices to iterate.
     * @param index Starting index in the indices vector.
     * @param strand Optional strand filter.
     */
    explicit Iterator(const IITree<int, GenomicFeature>* tree, std::span<const size_t> indices,
                      size_t index, const std::optional<GenomicStrand>& strand);

    [[nodiscard]] auto operator*() const -> reference;
    [[nodiscard]] auto operator->() const -> pointer;

    auto operator++() -> Iterator&;
    auto operator++(int) -> Iterator;

    friend auto operator==(const Iterator& lhs, const Iterator& rhs) -> bool {
        return lhs.current_index == rhs.current_index && lhs.indices.data() == rhs.indices.data();
    }

    friend auto operator!=(const Iterator& lhs, const Iterator& rhs) -> bool {
        return !(lhs == rhs);
    }

   private:
    const IITree<int, GenomicFeature>* tree;
    std::span<const size_t> indices;
    size_t current_index;
    std::optional<GenomicStrand> strand;
};

struct FeatureAnnotator::MergeInsertResult {
    std::string featureID;
    std::vector<std::string> mergedFeatureIDs;
};

}  // namespace annotation
