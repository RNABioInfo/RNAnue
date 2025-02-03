#pragma once

// Boost
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard
#include <optional>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

// Internal
#include "GenomicRegion.hpp"

namespace dataTypes {

struct GenomicFeature {
    std::string type;
    GenomicRegion genomicRegion;
    std::string featureID;
    std::optional<std::string> groupID;
    std::optional<std::string> geneName;
};

// TODO Rewrite as constructor of GenomicFeature
[[nodiscard]] inline auto asGenomicFeature(const GenomicRegion& region,
                                           std::string&& featureType = "supplementary_feature")
    -> GenomicFeature {
    namespace uuids = boost::uuids;
    const std::string uuid = uuids::to_string(uuids::random_generator()());

    return {.type = featureType,
            .genomicRegion = region,
            .featureID = uuid,
            .groupID = std::nullopt,
            .geneName = std::nullopt};
}

inline auto operator<<(std::ostream& ostream, const GenomicFeature& genomicFeature)
    -> std::ostream& {
    ostream << "type: " << genomicFeature.type << ", id: " << genomicFeature.featureID;

    if (genomicFeature.groupID) {
        ostream << ", groupID: " << genomicFeature.groupID.value();
    }

    if (genomicFeature.geneName) {
        ostream << ", geneName: " << genomicFeature.geneName.value();
    }

    ostream << ", region: " << genomicFeature.genomicRegion;

    return ostream;
}

using FeatureMap = std::unordered_map<int, std::vector<GenomicFeature>>;

}  // namespace dataTypes
