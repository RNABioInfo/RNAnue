#pragma once

// Internal
#include <string>
#include <utility>

#include "InteractionCluster.hpp"

namespace pipelines::analyze {

class PartiallyAnnotatedInteractionCluster : public InteractionCluster {
   public:
    PartiallyAnnotatedInteractionCluster(InteractionCluster&& cluster,
                                         std::optional<std::string> firstFeatureID,
                                         std::optional<std::string> secondFeatureID)
        : InteractionCluster(std::move(cluster)),
          firstFeatureID(std::move(firstFeatureID)),
          secondFeatureID(std::move(secondFeatureID)) {}

    [[nodiscard]] auto getFirstFeatureID() const -> const std::optional<std::string>& {
        return firstFeatureID;
    }

    [[nodiscard]] auto hasFirstFeatureID() const -> bool { return firstFeatureID.has_value(); }

    [[nodiscard]] auto hasSecondFeatureID() const -> bool { return secondFeatureID.has_value(); }

    [[nodiscard]] auto getSecondFeatureID() const -> const std::optional<std::string>& {
        return secondFeatureID;
    }

   private:
    std::optional<std::string> firstFeatureID;
    std::optional<std::string> secondFeatureID;
};

}  // namespace pipelines::analyze
