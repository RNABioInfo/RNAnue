#include "FeatureAnnotator.hpp"

// Standard
#include <algorithm>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iterator>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

// Internal
#include "FeatureParser.hpp"
#include "GenomicFeature.hpp"
#include "GenomicFeatureTreeMerger.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "IITree.hpp"
#include "SamRecord.hpp"

namespace annotation {
using std::unordered_set;

FeatureAnnotator::FeatureAnnotator(fs::path &featureFilePath,
                                   const ReferenceIDToIndexMap &referenceIDToIndex,
                                   const std::unordered_set<std::string> &includedFeatures,
                                   const std::string &featureIDFlag)
    : featureTreeMap(buildFeatureTreeMap(featureFilePath, referenceIDToIndex, includedFeatures,
                                         featureIDFlag)) {}

FeatureAnnotator::FeatureAnnotator(fs::path &featureFilePath,
                                   const ReferenceIDToIndexMap &referenceIDToIndex,
                                   const std::unordered_set<std::string> &includedFeatures)
    : featureTreeMap(buildFeatureTreeMap(featureFilePath, referenceIDToIndex, includedFeatures,
                                         std::nullopt)) {}

FeatureAnnotator::FeatureAnnotator(const FeatureMap &featureMap)
    : featureTreeMap(buildFeatureTreeMap(featureMap)) {}

auto FeatureAnnotator::buildFeatureTreeMap(const FeatureMap &featureMap) -> FeatureTreeMap {
    FeatureTreeMap newFeatureTreeMap;
    newFeatureTreeMap.reserve(featureMap.size());

    for (const auto &[referenceIDIndex, features] : featureMap) {
        IITree<int, dataTypes::GenomicFeature> tree;
        for (const auto &feature : features) {
            tree.add(feature.genomicRegion.getStart(), feature.genomicRegion.getEnd(), feature);
        }

        tree.index();

        newFeatureTreeMap.emplace(referenceIDIndex, std::move(tree));
    }

    return newFeatureTreeMap;
}

auto FeatureAnnotator::buildFeatureTreeMap(const fs::path &featureFilePath,
                                           const ReferenceIDToIndexMap &referenceIDToIndex,
                                           const std::unordered_set<std::string> &includedFeatures,
                                           const std::optional<std::string> &featureIDFlag)
    -> FeatureTreeMap {
    return buildFeatureTreeMap(
        FeatureParser(includedFeatures, featureIDFlag).parse(featureFilePath, referenceIDToIndex));
}

auto FeatureAnnotator::featureCount() const -> size_t {
    size_t count = 0;
    for (const auto &[_, tree] : featureTreeMap) {
        count += tree.size();
    }
    return count;
}

auto FeatureAnnotator::insert(const GenomicRegion &region) -> std::string {
    namespace uuids = boost::uuids;
    const std::string uuid = uuids::to_string(uuids::random_generator()());

    auto &tree = featureTreeMap[region.getReferenceIDIndex()];

    GenomicFeature feature{.type = "supplementary_feature",
                           .genomicRegion = region,
                           .featureID = uuid,
                           .groupID = std::nullopt,
                           .geneName = std::nullopt};

    tree.add(region.getStart(), region.getEnd(), feature);

    return uuid;
}

auto FeatureAnnotator::insertIndex(const GenomicRegion &region) -> std::string {
    namespace uuids = boost::uuids;
    const std::string uuid = uuids::to_string(uuids::random_generator()());

    auto &tree = featureTreeMap[region.getReferenceIDIndex()];

    GenomicFeature feature{.type = "supplementary_feature",
                           .genomicRegion = region,
                           .featureID = uuid,
                           .groupID = std::nullopt,
                           .geneName = std::nullopt};

    tree.add(region.getStart(), region.getEnd(), feature);

    tree.index();

    return uuid;
}

auto FeatureAnnotator::getOverlappingFeatures(const GenomicRegion &region,
                                              const GenomicOrientation orientation) const
    -> std::vector<GenomicFeature> {
    std::vector<GenomicFeature> features;
    auto iterator = featureTreeMap.find(region.getReferenceIDIndex());

    if (iterator != featureTreeMap.end()) {
        std::vector<size_t> indices;
        iterator->second.overlap(region.getStart(), region.getEnd(), indices);

        features.reserve(indices.size());

        for (const auto &index : indices) {
            const auto &feature = iterator->second.getData(index);

            if ((orientation == GenomicOrientation::BOTH) ||
                (orientation == GenomicOrientation::OPPOSITE &&
                 feature.genomicRegion.getStrand() == !region.getStrand()) ||
                (orientation == GenomicOrientation::SAME &&
                 feature.genomicRegion.getStrand() == region.getStrand())) {
                features.push_back(feature);
            }
        }
    }

    return features;
}

auto FeatureAnnotator::overlappingFeatureIterator(const GenomicRegion &region,
                                                  const GenomicOrientation orientation) const
    -> FeatureAnnotator::Results {
    std::vector<size_t> indices;

    auto iterator = featureTreeMap.find(region.getReferenceIDIndex());

    if (iterator == featureTreeMap.end()) [[unlikely]] {
        return {&iterator->second, indices, std::nullopt};
    }

    iterator->second.overlap(region.getStart(), region.getEnd(), indices);

    std::optional<dataTypes::GenomicStrand> strand = std::nullopt;

    if (orientation == GenomicOrientation::SAME) {
        strand = region.getStrand();
    } else if (orientation == GenomicOrientation::OPPOSITE) {
        strand = !region.getStrand();
    }

    return {&iterator->second, indices, strand};
}

auto FeatureAnnotator::getBestOverlappingFeature(const dataTypes::GenomicRegion &region,
                                                 const GenomicOrientation orientation) const
    -> std::optional<dataTypes::GenomicFeature> {
    auto overlapSizeWithRegion = [region](const dataTypes::GenomicFeature &feature) -> size_t {
        return std::min(region.getEnd(), feature.genomicRegion.getEnd()) -
               std::max(region.getStart(), feature.genomicRegion.getStart());
    };

    auto featureIterator = overlappingFeatureIterator(region, orientation);

    auto maxOverlapSizeElement =
        std::max_element(featureIterator.begin(), featureIterator.end(),
                         [&overlapSizeWithRegion](const dataTypes::GenomicFeature &lhs,
                                                  const dataTypes::GenomicFeature &rhs) {
                             return std::invoke(overlapSizeWithRegion, lhs) <
                                    std::invoke(overlapSizeWithRegion, rhs);
                         });

    if (maxOverlapSizeElement != featureIterator.end()) {
        return *maxOverlapSizeElement;
    }
    return std::nullopt;
}

auto FeatureAnnotator::getFeatureTreeMap() const -> const FeatureTreeMap & {
    return featureTreeMap;
}

auto FeatureAnnotator::mergeFeatures(const GenomicRegion &region, const int mergingTolerance)
    -> std::unordered_set<size_t> {
    auto featureTreeIterator = featureTreeMap.find(region.getReferenceIDIndex());

    if (featureTreeIterator == featureTreeMap.end()) [[unlikely]] {
        return {};
    }

    auto &featureTree = featureTreeIterator->second;

    // Merge blunt ended regions (no overlap)
    const int startSearchPosition = region.getStart() + mergingTolerance;
    const int endSearchPosition = region.getEnd() - mergingTolerance;

    std::vector<size_t> indices;
    featureTree.overlap(startSearchPosition, endSearchPosition, indices);

    std::erase_if(indices, [&region, &featureTree](size_t index) {
        return featureTree.getData(index).genomicRegion.getStrand() != region.getStrand();
    });

    if (indices.size() < 2) {
        return {};
    }

    // Find max end position.
    // Min position is already set to the start position of the first feature.
    auto maxEndIndex = std::ranges::max_element(indices, [&featureTree](size_t lhs, size_t rhs) {
        return featureTree.getData(lhs).genomicRegion.getEnd() <
               featureTree.getData(rhs).genomicRegion.getEnd();
    });

    int mergedStartPosition = featureTree.getData(indices[0]).genomicRegion.getStart();
    int mergedEndPosition = featureTree.getData(*maxEndIndex).genomicRegion.getEnd();

    if (mergedEndPosition > endSearchPosition) {
        auto newSearchRegion = GenomicRegion{region.getReferenceIDIndex(),
                                             {
                                                 .startPosition = mergedStartPosition,
                                                 .endPosition = mergedEndPosition,
                                             },
                                             region.getStrand()};

        return mergeFeatures(newSearchRegion, mergingTolerance);
    }

    featureTree.setIntervalEnd(indices[0], mergedEndPosition);
    featureTree.getData(indices[0]).genomicRegion.setEnd(mergedEndPosition);

    return {indices.begin() + 1, indices.end()};
}

void FeatureAnnotator::mergeIndexAllOverlappingFeatures(
    FeatureMergingParameters mergingParameters) {
    GenomicFeatureTreeMerger merger{mergingParameters};

    for (auto &tree : featureTreeMap) {
        tree.second.index();

        merger.merge(tree.second);
    }
}

auto FeatureAnnotator::getBestOverlappingFeature(const SamRecord &record,
                                                 const GenomicOrientation orientation) const
    -> std::optional<dataTypes::GenomicFeature> {
    auto region = GenomicRegion::fromSamRecord(record);

    if (!region.has_value()) {
        return std::nullopt;
    }

    return getBestOverlappingFeature(region.value(), orientation);
}

// Results and Iterator implementation
FeatureAnnotator::Results::Results(const IITree<int, dataTypes::GenomicFeature> *tree,
                                   const std::vector<size_t> &indices,
                                   const std::optional<dataTypes::GenomicStrand> strand)
    : tree(tree), indices(indices), strand(strand) {}

[[nodiscard]] auto FeatureAnnotator::Results::begin() const -> FeatureAnnotator::Results::Iterator {
    size_t startIndex = 0;
    if (strand.has_value()) {
        // Find the first index with the specified strand
        auto iterator = std::ranges::find_if(indices, [&](size_t index) {
            return tree->getData(index).genomicRegion.getStrand() == *strand;
        });
        startIndex =
            iterator != indices.end() ? std::distance(indices.begin(), iterator) : indices.size();
    }
    return Iterator(tree, indices, startIndex, strand);
}

[[nodiscard]] auto FeatureAnnotator::Results::end() const -> FeatureAnnotator::Results::Iterator {
    return Iterator(tree, indices, indices.size(), strand);
}

FeatureAnnotator::Results::Iterator::Iterator(const IITree<int, dataTypes::GenomicFeature> *tree,
                                              const std::vector<size_t> &indices, size_t index,
                                              const std::optional<dataTypes::GenomicStrand> &strand)
    : tree(tree), indices(indices), current_index(index), strand(strand) {}

[[nodiscard]] auto FeatureAnnotator::Results::Iterator::operator*() const
    -> FeatureAnnotator::Results::Iterator::reference {
    if (current_index >= indices.size()) {
        throw std::out_of_range("Iterator out of range");
    }
    return tree->getData(indices[current_index]);
}

[[nodiscard]] auto FeatureAnnotator::Results::Iterator::operator->() const
    -> FeatureAnnotator::Results::Iterator::pointer {
    if (current_index >= indices.size()) {
        throw std::out_of_range("Iterator out of range");
    }
    return &tree->getData(indices[current_index]);
}

auto FeatureAnnotator::Results::Iterator::operator++() -> FeatureAnnotator::Results::Iterator & {
    auto matchStrand = [&](size_t index) {
        if (strand.has_value()) {
            return tree->getData(indices[index]).genomicRegion.getStrand() == *strand;
        }
        return true;
    };

    ++current_index;

    while (current_index < indices.size() && !matchStrand(current_index)) {
        ++current_index;
    }

    return *this;
}

auto FeatureAnnotator::Results::Iterator::operator++(int) -> FeatureAnnotator::Results::Iterator {
    Iterator tmp = *this;
    ++(*this);
    return tmp;
}

}  // namespace annotation
