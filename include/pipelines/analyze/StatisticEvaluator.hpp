#pragma once

// Standard
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

// Boost
#include <boost/math/distributions/binomial.hpp>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "EvaluatedInteractionCluster.hpp"
#include "TranscriptContributionsByID.hpp"

namespace math = boost::math;

namespace pipelines::analyze {

class StatisticEvaluator {
   public:
    StatisticEvaluator() = delete;

    static auto evaluate(std::vector<AnnotatedInteractionCluster> &clusters,
                         const TranscriptContributionsByID &backgroundContributions,
                         double padjThreshold)
        -> std::vector<EvaluatedInteractionCluster>;

   private:
    double padjThreshold;

    static auto evaluatePValues(std::vector<AnnotatedInteractionCluster> &clusters,
                                const TranscriptContributionsByID &backgroundContributions)
        -> std::vector<EvaluatedInteractionCluster>;

    static auto evaluatePAdjValues(std::vector<EvaluatedInteractionCluster> &clusters,
                                   double pAdjThreshold)
        -> std::vector<EvaluatedInteractionCluster>;
};

}  // namespace pipelines::analyze
