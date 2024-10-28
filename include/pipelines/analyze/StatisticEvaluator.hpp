#pragma once

// Standard
#include <vector>

// Boost
#include <boost/math/distributions/binomial.hpp>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "EvaluatedInteractionCluster.hpp"

namespace math = boost::math;

namespace pipelines::analyze {

class StatisticEvaluator {
   public:
    StatisticEvaluator() = delete;

    static auto evaluate(std::vector<AnnotatedInteractionCluster> &clusters,
                         const std::unordered_map<std::string, size_t> &transcriptCounts,
                         size_t totalTranscriptCount, double padjThreshold)
        -> std::vector<EvaluatedInteractionCluster>;

   private:
    double padjThreshold;

    static auto getTranscriptProbabilities(
        const std::unordered_map<std::string, size_t> &transcriptCounts,
        size_t totalTranscriptCount) -> std::unordered_map<std::string, double>;

    static auto evaluatePValues(std::vector<AnnotatedInteractionCluster> &clusters,
                                const std::unordered_map<std::string, size_t> &transcriptCounts,
                                size_t totalTranscriptCount)
        -> std::vector<EvaluatedInteractionCluster>;

    static auto evaluatePAdjValues(std::vector<EvaluatedInteractionCluster> &clusters,
                                   double pAdjThreshold)
        -> std::vector<EvaluatedInteractionCluster>;
};

}  // namespace pipelines::analyze
