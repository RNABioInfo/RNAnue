#pragma once

// Standard
#include <cstddef>
#include <variant>

// Internal
#include "GenomicOrientation.hpp"
#include "GenomicStrandSpecificity.hpp"

namespace pipelines::analyze {

using namespace dataTypes;

struct ClusterOverlapToleranceMergeParameter {
    int tolerance{};

    explicit ClusterOverlapToleranceMergeParameter(int tolerance) noexcept : tolerance(tolerance) {}
};

[[nodiscard]] inline constexpr auto clusterDistanceToOverlapTolerance(int clusterDistance) noexcept
    -> int {
    // Distance 0 means blunt-ended intervals are allowed to merge, which is tolerance 1.
    return clusterDistance + 1;
}

struct ShortestClusterOverlapFractionMergeParameter {
    float overlapFraction{};

    explicit ShortestClusterOverlapFractionMergeParameter(float overlapFraction) noexcept
        : overlapFraction(overlapFraction) {}
};

using ClusteringMergeParameterVariant = std::variant<ClusterOverlapToleranceMergeParameter,
                                                     ShortestClusterOverlapFractionMergeParameter>;

struct ClusteringParameters {
    ClusteringMergeParameterVariant clusterMergeParameter;
    GenomicStrandSpecificity clusterMergingStrandSpecificity;
    double maxClusterSelfOverlapFraction;
    float minimumClusterTrascriptContribution;
    GenomicOrientation featureOrientation;
};

}  // namespace pipelines::analyze
