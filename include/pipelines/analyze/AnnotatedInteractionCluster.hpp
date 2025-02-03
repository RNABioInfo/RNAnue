#pragma once

// Standard
#include <cassert>
#include <optional>
#include <string>
#include <utility>

// Internal
#include "InteractionCluster.hpp"
#include "PartiallyAnnotatedInteractionCluster.hpp"

namespace pipelines::analyze {

class AnnotatedInteractionCluster : public InteractionCluster {
   public:
    AnnotatedInteractionCluster(const InteractionCluster& cluster, std::string firstFeatureID,
                                std::string secondFeatureID)
        : InteractionCluster(cluster),
          firstFeatureID(std::move(firstFeatureID)),
          secondFeatureID(std::move(secondFeatureID)) {}

    AnnotatedInteractionCluster(InteractionCluster&& cluster, std::string firstFeatureID,
                                std::string secondFeatureID)
        : InteractionCluster(std::move(cluster)),
          firstFeatureID(std::move(firstFeatureID)),
          secondFeatureID(std::move(secondFeatureID)) {}

    static auto fromPartiallyAnnotatedInteractionCluster(
        PartiallyAnnotatedInteractionCluster&& cluster, std::optional<std::string> firstFeatureID,
        std::optional<std::string> secondFeatureID) -> AnnotatedInteractionCluster {
        assert(cluster.hasFirstFeatureID() || firstFeatureID.has_value());
        assert(cluster.hasSecondFeatureID() || secondFeatureID.has_value());

        return {std::move(cluster), cluster.getFirstFeatureID().value_or(firstFeatureID.value()),
                cluster.getSecondFeatureID().value_or(secondFeatureID.value())};
    }

    [[nodiscard]] auto getFirstFeatureID() const -> const std::string& { return firstFeatureID; }

    [[nodiscard]] auto getSecondFeatureID() const -> const std::string& { return secondFeatureID; }

   private:
    std::string firstFeatureID;
    std::string secondFeatureID;
};

}  // namespace pipelines::analyze
