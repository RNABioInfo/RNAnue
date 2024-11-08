#pragma once

// Internal
#include "Orientation.hpp"

namespace pipelines::analyze {

struct ClusteringParameters {
    annotation::Orientation featureOrientation;
    double maxOverlapFraction;
    size_t minReadCount;
    int graceDistance;
};

}  // namespace pipelines::analyze
