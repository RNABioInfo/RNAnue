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
    size_t minimumClusterReadCount;
    GenomicOrientation featureOrientation;
};

}  // namespace pipelines::analyze
