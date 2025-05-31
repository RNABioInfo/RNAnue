#include "FeatureAnnotator.hpp"

// Standard
#include <algorithm>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <numeric>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/io/sam_file/input.hpp>

// Internal
#include "FeatureParser.hpp"
#include "GenomicFeature.hpp"
#include "GenomicFeatureTreeMerger.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "IITree.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "SamRecord.hpp"

namespace annotation {
using std::unordered_set;

auto loadReferenceIDToIndexMap(const std::vector<fs::path> &samFilePaths) -> ReferenceIDToIndexMap {
    ReferenceIDToIndexMap referenceToIndex;
    bool isFirst = true;

    for (const fs::path &inputPath : samFilePaths) {
        seqan3::sam_file_input splitsIn{inputPath, SamFieldIDs{}};

        auto &header = splitsIn.header();
        const std::deque<std::string> &referenceIDs = header.ref_ids();

        ReferenceIDToIndexMap sampleReferenceToIndex;
        sampleReferenceToIndex.reserve(referenceIDs.size());

        int index = 0;

        for (const auto &referenceID : referenceIDs) {
            sampleReferenceToIndex.emplace(referenceID, index++);
        }

        if (isFirst) {
            referenceToIndex = std::move(sampleReferenceToIndex);
            isFirst = false;
        } else if (referenceToIndex != sampleReferenceToIndex) {
            Logger::log<SourceLocation{}, LogLevel::ERROR>(
                "Samples seem to be mapped against different reference genomes, "
                "which is not "
                "supported.");
        }
    }

    return referenceToIndex;
}

FeatureAnnotator::FeatureAnnotator(fs::path &featureFilePath,
                                   const ReferenceIDToIndexMap &referenceIDToIndex,
                                   const std::unordered_set<std::string> &includedFeatures,
                                   const std::string &featureIDFlag) noexcept
    : featureTreeMap(buildFeatureTreeMap(featureFilePath, referenceIDToIndex, includedFeatures,
                                         featureIDFlag)) {}

FeatureAnnotator::FeatureAnnotator(fs::path &featureFilePath,
                                   const ReferenceIDToIndexMap &referenceIDToIndex,
                                   const std::unordered_set<std::string> &includedFeatures) noexcept
    : featureTreeMap(buildFeatureTreeMap(featureFilePath, referenceIDToIndex, includedFeatures,
                                         std::nullopt)) {}

FeatureAnnotator::FeatureAnnotator(const FeatureMap &featureMap) noexcept
    : featureTreeMap(buildFeatureTreeMap(featureMap)) {}

auto FeatureAnnotator::buildFeatureTreeMap(const FeatureMap &featureMap) -> FeatureTreeMap {
    FeatureTreeMap newFeatureTreeMap;
    newFeatureTreeMap.reserve(featureMap.size());

    for (const auto &[referenceIDIndex, features] : featureMap) {
        IITree<int, dataTypes::GenomicFeature> tree;
        for (const auto &feature : features) {
            tree.add(feature.getGenomicRegion().getStart(), feature.getGenomicRegion().getEnd(),
                     feature);
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

auto FeatureAnnotator::featureCount() const noexcept -> size_t {
    return std::accumulate(
        featureTreeMap.begin(), featureTreeMap.end(), size_t{0},
        [](size_t total, const auto &pair) noexcept { return total + pair.second.size(); });
}

auto FeatureAnnotator::insertIndex(const GenomicRegion &region) -> std::string {
    namespace uuids = boost::uuids;
    const std::string uuid = uuids::to_string(uuids::random_generator()());

    auto &tree = featureTreeMap[region.getReferenceIDIndex()];

    GenomicFeature feature{"supplementary_feature", region, uuid, std::nullopt, std::nullopt};

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
                 feature.getGenomicRegion().getStrand() == !region.getStrand()) ||
                (orientation == GenomicOrientation::SAME &&
                 feature.getGenomicRegion().getStrand() == region.getStrand())) {
                features.push_back(feature);
            }
        }
    }

    return features;
}

auto FeatureAnnotator::overlappingFeatureIt(const GenomicRegion &region,
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

auto FeatureAnnotator::getBestOverlappingFeature(const GenomicRegion &region,
                                                 const GenomicOrientation orientation) const
    -> std::optional<GenomicFeature> {
    auto featureIterator = overlappingFeatureIt(region, orientation);

    auto maxOverlapSizeElement = std::max_element(
        featureIterator.begin(), featureIterator.end(),
        [&region](const GenomicFeature &lhs, const GenomicFeature &rhs) {
            return lhs.getGenomicRegion().overlap(region) < rhs.getGenomicRegion().overlap(region);
        });

    if (maxOverlapSizeElement != featureIterator.end()) {
        return *maxOverlapSizeElement;
    }
    return std::nullopt;
}

auto FeatureAnnotator::getBestOverlappingFeature(const SamRecord &record,
                                                 const GenomicOrientation orientation) const
    -> std::optional<GenomicFeature> {
    auto region = GenomicRegion::fromSamRecord(record);

    if (!region.has_value()) {
        return std::nullopt;
    }

    return getBestOverlappingFeature(region.value(), orientation);
}

auto FeatureAnnotator::getFeatureTreeMap() const -> const FeatureTreeMap & {
    return featureTreeMap;
}

void FeatureAnnotator::mergeAllOverlappingFeatures(FeatureMergingParameters mergingParameters) {
    GenomicFeatureTreeMerger merger{mergingParameters};

    for (auto &[refIndex, tree] : featureTreeMap) {
        tree.index();
        merger.merge(tree);
    }
}

void FeatureAnnotator::debugOutputAllFeatures() const {
    for (const auto &[featureID, tree] : featureTreeMap) {
        std::cout << "Chromosome ID: " << featureID << "\n";
        for (const auto &feature : tree.intervals()) {
            std::cout << "Start: " << feature.start << ", End: " << feature.end
                      << "; Max:" << feature.max << "; ";
            std::cout << feature.data << "\n";
        }
    }
};

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
            return tree->getData(index).getGenomicRegion().getStrand() == *strand;
        });
        startIndex =
            iterator != indices.end() ? std::distance(indices.begin(), iterator) : indices.size();
    }
    return Iterator(tree, std::span<const size_t>(indices), startIndex, strand);
}

[[nodiscard]] auto FeatureAnnotator::Results::end() const -> FeatureAnnotator::Results::Iterator {
    return Iterator(tree, indices, indices.size(), strand);
}

FeatureAnnotator::Results::Iterator::Iterator(const IITree<int, dataTypes::GenomicFeature> *tree,
                                              std::span<const size_t> indices, size_t index,
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
            return tree->getData(indices[index]).getGenomicRegion().getStrand() == *strand;
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
