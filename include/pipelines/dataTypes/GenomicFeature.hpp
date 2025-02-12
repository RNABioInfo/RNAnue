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
#include <utility>
#include <vector>

// Internal
#include "GenomicRegion.hpp"

namespace dataTypes {

struct GenomicFeature {
   public:
    GenomicFeature(std::string type, GenomicRegion genomicRegion, std::string featureID,
                   std::optional<std::string> groupID, std::optional<std::string> geneName)
        : type(std::move(type)),
          genomicRegion(genomicRegion),
          featureID(std::move(featureID)),
          groupID(std::move(groupID)),
          geneName(std::move(geneName)) {}

    // Getters
    [[nodiscard]] auto getType() const noexcept -> const std::string& { return type; }
    [[nodiscard]] auto getGenomicRegion() const noexcept -> GenomicRegion { return genomicRegion; }
    [[nodiscard]] auto getGenomicRegion() noexcept -> GenomicRegion& { return genomicRegion; }
    [[nodiscard]] auto getID() const noexcept -> const std::string& { return featureID; }
    [[nodiscard]] auto getGroupID() const noexcept -> const std::optional<std::string>& {
        return groupID;
    }
    [[nodiscard]] auto getGeneName() const noexcept -> const std::optional<std::string>& {
        return geneName;
    }

    /**
     * @brief Returns the annotation identifier for this genomic feature.
     *
     * If a group ID is provided, it is used as the annotation ID; otherwise, the feature ID is
     * returned.
     *
     * @return std::string The annotation identifier.
     */
    [[nodiscard]] auto getAnnotationID() const noexcept -> std::string {
        return groupID.value_or(featureID);
    }

   private:
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

    return {featureType, region, uuid, std::nullopt, std::nullopt};
}

inline auto operator<<(std::ostream& ostream, const GenomicFeature& genomicFeature)
    -> std::ostream& {
    ostream << "type: " << genomicFeature.getType() << ", id: " << genomicFeature.getID();

    if (genomicFeature.getGroupID()) {
        ostream << ", groupID: " << genomicFeature.getGroupID().value();
    }

    if (genomicFeature.getGeneName()) {
        ostream << ", geneName: " << genomicFeature.getGeneName().value();
    }

    ostream << ", region: " << genomicFeature.getGenomicRegion();

    return ostream;
}

using FeatureMap = std::unordered_map<int, std::vector<GenomicFeature>>;

}  // namespace dataTypes
