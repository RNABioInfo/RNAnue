#pragma once

// Standard
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <variant>
#include <vector>

// Internal
#include "ClusteringParameters.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "IITree.hpp"
#include "InteractionCluster.hpp"
#include "VariantOverload.hpp"

namespace pipelines::analyze {

class InteractionClusterComponentBuilder {
   public:
    [[nodiscard]] static auto groupClusters(std::vector<InteractionCluster>&& clusters,
                                            GenomicStrandSpecificity strandSpecificity)
        -> std::vector<std::vector<InteractionCluster>> {
        if (clusters.empty()) {
            return {};
        }

        std::ranges::sort(clusters, [strandSpecificity](const InteractionCluster& lhs,
                                                        const InteractionCluster& rhs) {
            return clusterLess(lhs, rhs, strandSpecificity);
        });

        std::vector<std::vector<InteractionCluster>> groups;
        std::vector<InteractionCluster> currentGroup;
        currentGroup.emplace_back(std::move(clusters.front()));

        for (size_t index = 1; index < clusters.size(); ++index) {
            if (!sameGroup(currentGroup.front(), clusters[index], strandSpecificity)) {
                groups.emplace_back(std::move(currentGroup));
                currentGroup = {};
            }

            currentGroup.emplace_back(std::move(clusters[index]));
        }

        groups.emplace_back(std::move(currentGroup));
        return groups;
    }

    [[nodiscard]] static auto mergeGroup(std::vector<InteractionCluster>&& clusters,
                                         const ClusteringParameters& parameters)
        -> std::vector<InteractionCluster> {
        if (clusters.size() < 2) {
            return std::move(clusters);
        }

        IITree<int32_t, size_t> firstArmIndex;
        IITree<int32_t, size_t> secondArmIndex;

        for (size_t index = 0; index < clusters.size(); ++index) {
            const auto& firstSegment = clusters[index].getFirstSegment();
            const auto& secondSegment = clusters[index].getSecondSegment();
            firstArmIndex.add(firstSegment.getStart(), firstSegment.getEnd(), index);
            secondArmIndex.add(secondSegment.getStart(), secondSegment.getEnd(), index);
        }

        firstArmIndex.index();
        secondArmIndex.index();

        UnionFind components{clusters.size()};
        std::vector<size_t> firstArmHits;
        std::vector<size_t> secondArmHits;

        for (size_t index = 0; index < clusters.size(); ++index) {
            const auto firstQuery = queryRange(clusters[index].getFirstSegment(), parameters);
            const auto secondQuery = queryRange(clusters[index].getSecondSegment(), parameters);

            if (!firstQuery.valid || !secondQuery.valid) {
                continue;
            }

            firstArmIndex.overlap(firstQuery.start, firstQuery.end, firstArmHits);
            secondArmIndex.overlap(secondQuery.start, secondQuery.end, secondArmHits);

            const bool useFirstArmHits = firstArmHits.size() <= secondArmHits.size();
            const auto& selectedHits = useFirstArmHits ? firstArmHits : secondArmHits;
            const auto& selectedIndex = useFirstArmHits ? firstArmIndex : secondArmIndex;

            for (const size_t treeHitIndex : selectedHits) {
                const size_t candidateIndex = selectedIndex.getData(treeHitIndex);

                if (candidateIndex >= index) {
                    continue;
                }

                if (clustersOverlap(clusters[candidateIndex], clusters[index], parameters)) {
                    components.unite(candidateIndex, index);
                }
            }
        }

        return materializeComponents(std::move(clusters), components);
    }

    [[nodiscard]] static auto clustersOverlap(const InteractionCluster& cluster1,
                                              const InteractionCluster& cluster2,
                                              const ClusteringParameters& parameters) noexcept
        -> bool {
        return std::visit(
            overloaded{[&](const ClusterOverlapToleranceMergeParameter& tolerance) {
                           return cluster1.overlapsWithTolerance(
                               cluster2, parameters.clusterMergingStrandSpecificity,
                               tolerance.tolerance);
                       },
                       [&](const ShortestClusterOverlapFractionMergeParameter& overlapFraction) {
                           return cluster1.overlapsWithShortestSegmentFraction(
                               cluster2, parameters.clusterMergingStrandSpecificity,
                               overlapFraction.overlapFraction);
                       }},
            parameters.clusterMergeParameter);
    }

   private:
    struct QueryRange {
        int32_t start{};
        int32_t end{};
        bool valid{false};
    };

    class UnionFind {
       public:
        explicit UnionFind(size_t size) : parents(size), ranks(size, 0) {
            std::iota(parents.begin(), parents.end(), size_t{0});
        }

        [[nodiscard]] auto find(size_t value) -> size_t {
            if (parents[value] != value) {
                parents[value] = find(parents[value]);
            }

            return parents[value];
        }

        void unite(size_t lhs, size_t rhs) {
            size_t lhsRoot = find(lhs);
            size_t rhsRoot = find(rhs);

            if (lhsRoot == rhsRoot) {
                return;
            }

            if (ranks[lhsRoot] < ranks[rhsRoot]) {
                std::swap(lhsRoot, rhsRoot);
            }

            parents[rhsRoot] = lhsRoot;

            if (ranks[lhsRoot] == ranks[rhsRoot]) {
                ++ranks[lhsRoot];
            }
        }

       private:
        std::vector<size_t> parents;
        std::vector<size_t> ranks;
    };

    [[nodiscard]] static auto queryRange(const GenomicRegion& region,
                                         const ClusteringParameters& parameters) noexcept
        -> QueryRange {
        return std::visit(
            overloaded{[&](const ClusterOverlapToleranceMergeParameter& tolerance) {
                           return toleranceQueryRange(region, tolerance.tolerance);
                       },
                       [&](const ShortestClusterOverlapFractionMergeParameter&) {
                           return exactIntervalQueryRange(region);
                       }},
            parameters.clusterMergeParameter);
    }

    [[nodiscard]] static auto toleranceQueryRange(const GenomicRegion& region,
                                                  int tolerance) noexcept -> QueryRange {
        const auto start = clampToInt32(static_cast<int64_t>(region.getStart()) - tolerance);
        const auto end = clampToInt32(static_cast<int64_t>(region.getEnd()) + tolerance);

        return {.start = start, .end = end, .valid = start < end};
    }

    [[nodiscard]] static auto exactIntervalQueryRange(const GenomicRegion& region) noexcept
        -> QueryRange {
        return {.start = region.getStart(), .end = region.getEnd(),
                .valid = region.getStart() < region.getEnd()};
    }

    [[nodiscard]] static auto clampToInt32(int64_t value) noexcept -> int32_t {
        if (value < std::numeric_limits<int32_t>::min()) {
            return std::numeric_limits<int32_t>::min();
        }

        if (value > std::numeric_limits<int32_t>::max()) {
            return std::numeric_limits<int32_t>::max();
        }

        return static_cast<int32_t>(value);
    }

    [[nodiscard]] static auto materializeComponents(std::vector<InteractionCluster>&& clusters,
                                                    UnionFind& components)
        -> std::vector<InteractionCluster> {
        constexpr size_t noComponent = std::numeric_limits<size_t>::max();

        std::vector<std::vector<size_t>> componentMembers;
        std::vector<size_t> componentIndexByRoot(clusters.size(), noComponent);

        for (size_t index = 0; index < clusters.size(); ++index) {
            const size_t root = components.find(index);

            if (componentIndexByRoot[root] == noComponent) {
                componentIndexByRoot[root] = componentMembers.size();
                componentMembers.emplace_back();
            }

            componentMembers[componentIndexByRoot[root]].push_back(index);
        }

        std::vector<InteractionCluster> mergedClusters;
        mergedClusters.reserve(componentMembers.size());

        for (const auto& members : componentMembers) {
            InteractionCluster mergedCluster{std::move(clusters[members.front()])};

            for (size_t memberIndex = 1; memberIndex < members.size(); ++memberIndex) {
                mergedCluster.absorbValidatedComponentMember(clusters[members[memberIndex]]);
            }

            mergedClusters.emplace_back(std::move(mergedCluster));
        }

        return mergedClusters;
    }

    [[nodiscard]] static auto sameGroup(const InteractionCluster& lhs,
                                        const InteractionCluster& rhs,
                                        GenomicStrandSpecificity strandSpecificity) noexcept
        -> bool {
        if (lhs.getFirstSegment().getReferenceIDIndex() !=
                rhs.getFirstSegment().getReferenceIDIndex() ||
            lhs.getSecondSegment().getReferenceIDIndex() !=
                rhs.getSecondSegment().getReferenceIDIndex()) {
            return false;
        }

        if (strandSpecificity == GenomicStrandSpecificity::UNSPECIFIC) {
            return true;
        }

        return lhs.getFirstSegment().getStrand() == rhs.getFirstSegment().getStrand() &&
               lhs.getSecondSegment().getStrand() == rhs.getSecondSegment().getStrand();
    }

    [[nodiscard]] static auto clusterLess(const InteractionCluster& lhs,
                                          const InteractionCluster& rhs,
                                          GenomicStrandSpecificity strandSpecificity) noexcept
        -> bool {
        if (const int groupComparison = compareGroupKey(lhs, rhs, strandSpecificity);
            groupComparison != 0) {
            return groupComparison < 0;
        }

        return compareClusterCoordinates(lhs, rhs) < 0;
    }

    [[nodiscard]] static auto compareGroupKey(const InteractionCluster& lhs,
                                              const InteractionCluster& rhs,
                                              GenomicStrandSpecificity strandSpecificity) noexcept
        -> int {
        if (const int comparison = compare(lhs.getFirstSegment().getReferenceIDIndex(),
                                           rhs.getFirstSegment().getReferenceIDIndex());
            comparison != 0) {
            return comparison;
        }

        if (const int comparison = compare(lhs.getSecondSegment().getReferenceIDIndex(),
                                           rhs.getSecondSegment().getReferenceIDIndex());
            comparison != 0) {
            return comparison;
        }

        if (strandSpecificity == GenomicStrandSpecificity::UNSPECIFIC) {
            return 0;
        }

        if (const int comparison =
                compareStrand(lhs.getFirstSegment().getStrand(), rhs.getFirstSegment().getStrand());
            comparison != 0) {
            return comparison;
        }

        return compareStrand(lhs.getSecondSegment().getStrand(),
                             rhs.getSecondSegment().getStrand());
    }

    [[nodiscard]] static auto compareClusterCoordinates(const InteractionCluster& lhs,
                                                        const InteractionCluster& rhs) noexcept
        -> int {
        if (const int comparison =
                compare(lhs.getFirstSegment().getStart(), rhs.getFirstSegment().getStart());
            comparison != 0) {
            return comparison;
        }

        if (const int comparison =
                compare(lhs.getFirstSegment().getEnd(), rhs.getFirstSegment().getEnd());
            comparison != 0) {
            return comparison;
        }

        if (const int comparison =
                compare(lhs.getSecondSegment().getStart(), rhs.getSecondSegment().getStart());
            comparison != 0) {
            return comparison;
        }

        if (const int comparison =
                compare(lhs.getSecondSegment().getEnd(), rhs.getSecondSegment().getEnd());
            comparison != 0) {
            return comparison;
        }

        if (const int comparison =
                compareStrand(lhs.getFirstSegment().getStrand(), rhs.getFirstSegment().getStrand());
            comparison != 0) {
            return comparison;
        }

        return compareStrand(lhs.getSecondSegment().getStrand(),
                             rhs.getSecondSegment().getStrand());
    }

    template <typename T>
    [[nodiscard]] static auto compare(const T& lhs, const T& rhs) noexcept -> int {
        if (lhs < rhs) {
            return -1;
        }

        if (rhs < lhs) {
            return 1;
        }

        return 0;
    }

    [[nodiscard]] static auto compareStrand(GenomicStrand lhs, GenomicStrand rhs) noexcept -> int {
        return compare(static_cast<char>(lhs), static_cast<char>(rhs));
    }
};

}  // namespace pipelines::analyze
