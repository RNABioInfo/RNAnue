#pragma once

// Boost
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard
#include <cstddef>
#include <deque>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Internal
#include "GenomicFeature.hpp"
#include "GenomicFeatureTreeMerger.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "IITree.hpp"
#include "Logger.hpp"
#include "SamRecord.hpp"
#include "seqan3/io/sam_file/input.hpp"

namespace annotation {

using namespace dataTypes;

using ReferenceIDIndex = int;
using ReferenceIDToIndexMap = std::unordered_map<std::string, int>;
using FeatureTreeMap = std::unordered_map<ReferenceIDIndex, IITree<int, GenomicFeature>>;

namespace fs = std::filesystem;

[[nodiscard]] inline auto loadReferenceIDToIndexMap(const std::vector<fs::path>& samFilePaths)
    -> ReferenceIDToIndexMap {
    ReferenceIDToIndexMap referenceToIndex;

    bool isFirst = true;

    for (const fs::path& inputPath : samFilePaths) {
        seqan3::sam_file_input splitsIn{inputPath, SamFieldIDs{}};

        auto& header = splitsIn.header();
        const std::deque<std::string>& referenceIDs = header.ref_ids();

        ReferenceIDToIndexMap sampleReferenceToIndex;
        sampleReferenceToIndex.reserve(referenceIDs.size());

        for (int index = 0; const auto& referenceID : referenceIDs) {
            sampleReferenceToIndex.emplace(referenceID, index++);
        }

        if (isFirst) {
            referenceToIndex = std::move(sampleReferenceToIndex);
            isFirst = false;
        } else if (referenceToIndex != sampleReferenceToIndex) {
            Logger::log<SourceLocation{}, LogLevel::ERROR>(
                "Samples seem to be mapped against different reference genomes, which is not "
                "supported.");
        }
    }

    return referenceToIndex;
}

class FeatureAnnotator {
   public:
    FeatureAnnotator(fs::path& featureFilePath, const ReferenceIDToIndexMap& referenceIDToIndex,
                     const std::unordered_set<std::string>& includedFeatures,
                     const std::string& featureIDFlag);
    FeatureAnnotator(fs::path& featureFilePath, const ReferenceIDToIndexMap& referenceIDToIndex,
                     const std::unordered_set<std::string>& includedFeatures);

    explicit FeatureAnnotator(const FeatureMap& featureMap);

    explicit FeatureAnnotator() = default;

    FeatureAnnotator(const FeatureAnnotator&) = default;
    FeatureAnnotator(FeatureAnnotator&&) = default;
    auto operator=(const FeatureAnnotator&) -> FeatureAnnotator& = default;
    auto operator=(FeatureAnnotator&&) -> FeatureAnnotator& = delete;
    ~FeatureAnnotator() = default;

    class Results;
    struct MergeInsertResult;

    [[nodiscard]] auto featureCount() const -> size_t;

    auto insertIndex(const GenomicRegion& region) -> std::string;
    auto insert(const GenomicRegion& region) -> std::string;

    [[nodiscard]] auto getOverlappingFeatures(const GenomicRegion& region,
                                              GenomicOrientation orientation) const
        -> std::vector<GenomicFeature>;
    [[nodiscard]] auto overlappingFeatureIterator(const GenomicRegion& region,
                                                  GenomicOrientation orientation) const -> Results;
    [[nodiscard]] auto getBestOverlappingFeature(const GenomicRegion& region,
                                                 GenomicOrientation orientation) const
        -> std::optional<GenomicFeature>;
    [[nodiscard]] auto getBestOverlappingFeature(const SamRecord& record,
                                                 GenomicOrientation orientation) const
        -> std::optional<GenomicFeature>;

    [[nodiscard]] auto getFeatureTreeMap() const -> const FeatureTreeMap&;

    void mergeIndexAllOverlappingFeatures(FeatureMergingParameters mergingParameters);

    void printAllFeatures() const {
        for (const auto& [featureID, tree] : featureTreeMap) {
            std::cout << "Chromosome ID: " << featureID << "\n";
            for (const auto& feature : tree.intervals()) {
                std::cout << "Start: " << feature.start << ", End: " << feature.end
                          << "; Max:" << feature.max << "; ";
                std::cout << feature.data << "\n";
            }
        }
    };

   private:
    FeatureTreeMap featureTreeMap;

    static auto buildFeatureTreeMap(const fs::path& featureFilePath,
                                    const ReferenceIDToIndexMap& referenceIDToIndex,
                                    const std::vector<std::string>& includedFeatures,
                                    const std::optional<std::string>& featureIDFlag)
        -> FeatureTreeMap;
    static auto buildFeatureTreeMap(const fs::path& featureFilePath,
                                    const ReferenceIDToIndexMap& referenceIDToIndex,
                                    const std::unordered_set<std::string>& includedFeatures,
                                    const std::optional<std::string>& featureIDFlag)
        -> FeatureTreeMap;
    static auto buildFeatureTreeMap(const FeatureMap& featureMap) -> FeatureTreeMap;

    [[nodiscard]] auto mergeFeatures(const GenomicRegion& region, int tolerance)
        -> std::unordered_set<size_t>;
};

class FeatureAnnotator::Results {
   public:
    Results(const IITree<int, GenomicFeature>* tree, const std::vector<size_t>& indices,
            std::optional<GenomicStrand> strand);

    Results() = delete;

    struct Iterator;

    [[nodiscard]] auto begin() const -> Iterator;
    [[nodiscard]] auto end() const -> Iterator;

   private:
    const IITree<int, GenomicFeature>* tree;
    std::vector<size_t> indices;
    std::optional<GenomicStrand> strand;
};

struct FeatureAnnotator::Results::Iterator {
    using value_type = GenomicFeature;
    using difference_type = std::ptrdiff_t;
    using reference = const GenomicFeature&;
    using pointer = const GenomicFeature*;
    using iterator_category = std::forward_iterator_tag;

    explicit Iterator(const IITree<int, GenomicFeature>* tree, const std::vector<size_t>& indices,
                      size_t index, const std::optional<GenomicStrand>& strand);

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
    const IITree<int, GenomicFeature>* tree;
    std::vector<size_t> indices;
    size_t current_index;
    std::optional<GenomicStrand> strand;
};

struct FeatureAnnotator::MergeInsertResult {
    std::string featureID;
    std::vector<std::string> mergedFeatureIDs;
};

}  // namespace annotation
