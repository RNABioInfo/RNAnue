#pragma once

// Boost
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard
#include <deque>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicRegion.hpp"
#include "IITree.hpp"
#include "Orientation.hpp"
#include "SamRecord.hpp"

using namespace dataTypes;

namespace annotation {

using FeatureTreeMap = std::unordered_map<std::string, IITree<int, dataTypes::GenomicFeature>>;
namespace fs = std::filesystem;

class FeatureAnnotator {
   public:
    FeatureAnnotator(fs::path& featureFilePath,
                     const std::unordered_set<std::string>& includedFeatures,
                     const std::string& featureIDFlag);
    FeatureAnnotator(fs::path& featureFilePath,
                     const std::unordered_set<std::string>& includedFeatures);

    explicit FeatureAnnotator(const dataTypes::FeatureMap& featureMap);

    explicit FeatureAnnotator() = default;
    FeatureAnnotator(const FeatureAnnotator&) = default;
    FeatureAnnotator(FeatureAnnotator&&) = default;
    auto operator=(const FeatureAnnotator&) -> FeatureAnnotator& = default;
    auto operator=(FeatureAnnotator&&) -> FeatureAnnotator& = delete;
    ~FeatureAnnotator() = default;

    class Results;
    struct MergeInsertResult;

    [[nodiscard]] auto featureCount() const -> size_t;

    auto insertIndex(const dataTypes::GenomicRegion& region) -> std::string;
    auto insert(const dataTypes::GenomicRegion& region) -> std::string;

    auto mergeInsertIndex(const dataTypes::GenomicRegion& region, int graceDistance)
        -> MergeInsertResult;

    auto getOverlappingFeatures(const dataTypes::GenomicRegion& region, Orientation orientation)
        -> std::vector<dataTypes::GenomicFeature>;
    [[nodiscard]] auto overlappingFeatureIterator(const dataTypes::GenomicRegion& region,
                                                  Orientation orientation) const -> Results;
    [[nodiscard]] auto getBestOverlappingFeature(const dataTypes::GenomicRegion& region,
                                                 Orientation orientation) const
        -> std::optional<dataTypes::GenomicFeature>;
    [[nodiscard]] auto getBestOverlappingFeature(const SamRecord& record,
                                                 const std::deque<std::string>& referenceIDs,
                                                 Orientation orientation) const
        -> std::optional<dataTypes::GenomicFeature>;

    [[nodiscard]] auto getFeatureTreeMap() const -> const FeatureTreeMap&;

    void mergeIndexAllOverlappingFeatures(int minOverlap);

    void printAllFeatures() const {
        for (const auto& [featureID, tree] : featureTreeMap) {
            std::cout << "Chromosome ID: " << featureID << std::endl;
            for (const auto& feature : tree.intervals()) {
                std::cout << feature.data.id << feature.data.referenceID << ":"
                          << feature.data.startPosition << "-" << feature.data.endPosition << "\n";
            }
        }
    };

   private:
    FeatureTreeMap featureTreeMap;

    static auto buildFeatureTreeMap(const fs::path& featureFilePath,
                                    const std::vector<std::string>& includedFeatures,
                                    const std::optional<std::string>& featureIDFlag)
        -> FeatureTreeMap;
    static auto buildFeatureTreeMap(const fs::path& featureFilePath,
                                    const std::unordered_set<std::string>& includedFeatures,
                                    const std::optional<std::string>& featureIDFlag)
        -> FeatureTreeMap;
    static auto buildFeatureTreeMap(const dataTypes::FeatureMap& featureMap) -> FeatureTreeMap;

    [[nodiscard]] auto mergeFeatures(const dataTypes::GenomicRegion& region, int minOverlap)
        -> std::unordered_set<size_t>;
};

class FeatureAnnotator::Results {
   public:
    Results(const IITree<int, dataTypes::GenomicFeature>* tree, const std::vector<size_t>& indices,
            std::optional<dataTypes::GenomicStrand> strand);

    Results() = delete;

    struct Iterator;

    [[nodiscard]] auto begin() const -> Iterator;
    [[nodiscard]] auto end() const -> Iterator;

   private:
    const IITree<int, dataTypes::GenomicFeature>* tree;
    std::vector<size_t> indices;
    std::optional<dataTypes::GenomicStrand> strand;
};

struct FeatureAnnotator::Results::Iterator {
    using value_type = dataTypes::GenomicFeature;
    using difference_type = std::ptrdiff_t;
    using reference = const dataTypes::GenomicFeature&;
    using pointer = const dataTypes::GenomicFeature*;
    using iterator_category = std::forward_iterator_tag;

    explicit Iterator(const IITree<int, dataTypes::GenomicFeature>* tree,
                      const std::vector<size_t>& indices, size_t index,
                      const std::optional<dataTypes::GenomicStrand>& strand);

    [[nodiscard]] auto operator*() const -> reference;
    [[nodiscard]] auto operator->() const -> pointer;

    auto operator++() -> Iterator&;
    auto operator++(int) -> Iterator;

    friend auto operator==(const Iterator& lhs, const Iterator& rhs) -> bool {
        return lhs.current_index == rhs.current_index;
    }

    friend auto operator!=(const Iterator& lhs, const Iterator& rhs) -> bool {
        return lhs.current_index != rhs.current_index;
    }

   private:
    const IITree<int, dataTypes::GenomicFeature>* tree;
    std::vector<size_t> indices;
    size_t current_index;
    std::optional<dataTypes::GenomicStrand> strand;
};

struct FeatureAnnotator::MergeInsertResult {
    std::string featureID;
    std::vector<std::string> mergedFeatureIDs;
};

}  // namespace annotation
