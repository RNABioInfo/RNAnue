#include "StatisticEvaluator.hpp"

// Standard
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Boost
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/detail/derived_accessors.hpp>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "EvaluatedInteractionCluster.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "TranscriptContributionsByID.hpp"

namespace pipelines::analyze {

auto StatisticEvaluator::evaluate(std::vector<AnnotatedInteractionCluster> &clusters,
                                  const TranscriptContributionsByID &transcriptCounts,
                                  float totalTranscriptContribution, double padjThreshold)
    -> std::vector<EvaluatedInteractionCluster> {
    auto evaluatedClusters =
        evaluatePValues(clusters, transcriptCounts, totalTranscriptContribution);
    return evaluatePAdjValues(evaluatedClusters, padjThreshold);
}

auto StatisticEvaluator::getTranscriptProbabilities(
    const TranscriptContributionsByID &transcriptCounts, const float totalTranscriptContribution)
    -> std::unordered_map<std::string, float> {
    std::unordered_map<std::string, float> transcriptProbabilities;
    transcriptProbabilities.reserve(transcriptCounts.size());

    for (const auto &[transcriptID, count] : transcriptCounts) {
        transcriptProbabilities[transcriptID] = count / totalTranscriptContribution;
    }

    return transcriptProbabilities;
}

auto StatisticEvaluator::evaluatePValues(std::vector<AnnotatedInteractionCluster> &clusters,
                                         const TranscriptContributionsByID &transcriptCounts,
                                         const float totalTranscriptContribution)
    -> std::vector<EvaluatedInteractionCluster> {
    const std::unordered_map<std::string, float> transcriptProbabilities =
        getTranscriptProbabilities(transcriptCounts, totalTranscriptContribution);

    std::vector<EvaluatedInteractionCluster> evaluatedClusters;
    evaluatedClusters.reserve(clusters.size());

    double combinedProbability = 0.0;

    std::vector<double> ligationByChanceProbabilities;
    ligationByChanceProbabilities.reserve(clusters.size());

    for (auto &cluster : clusters) {
        const auto &firstTranscriptID = cluster.getFirstFeatureID();
        const auto &secondTranscriptID = cluster.getSecondFeatureID();

        auto firstIt = transcriptProbabilities.find(firstTranscriptID);
        auto secondIt = transcriptProbabilities.find(secondTranscriptID);

        if (firstIt == transcriptProbabilities.end() || secondIt == transcriptProbabilities.end()) {
            Logger::log<LogLevel::WARNING>(
                "Could not find transcript probabilities for cluster with transcripts: ",
                firstTranscriptID, ", ", secondTranscriptID);
            continue;
        }

        // Compute the ligation probability
        const double firstTranscriptProbability = firstIt->second;
        const double secondTranscriptProbability = secondIt->second;

        const double ligationByChanceProbability =
            (firstTranscriptID == secondTranscriptID)
                ? firstTranscriptProbability * secondTranscriptProbability
                : 2 * firstTranscriptProbability * secondTranscriptProbability;

        combinedProbability += ligationByChanceProbability;

        ligationByChanceProbabilities.emplace_back(ligationByChanceProbability);
    }

    for (size_t i = 0; i < clusters.size(); ++i) {
        auto &cluster = clusters[i];

        double normalizedLigationByChanceProbability =
            ligationByChanceProbabilities[i] / combinedProbability;

        if (std::isnan(normalizedLigationByChanceProbability) ||
            normalizedLigationByChanceProbability < 0.0 ||
            normalizedLigationByChanceProbability > 1.0) {
            Logger::log<LogLevel::DEBUG>("Skipping record due ligation by chance prob invalid.");
        }

        const auto binomialDistribution =
            math::binomial_distribution(static_cast<double>(totalTranscriptContribution),
                                        normalizedLigationByChanceProbability);

        const double pValue =
            1 - math::cdf(binomialDistribution, cluster.getTranscriptContribution() - 1);

        evaluatedClusters.emplace_back(std::move(cluster), pValue);
    }

    return evaluatedClusters;
}

auto StatisticEvaluator::evaluatePAdjValues(std::vector<EvaluatedInteractionCluster> &clusters,
                                            const double pAdjThreshold)
    -> std::vector<EvaluatedInteractionCluster> {
    std::vector<std::pair<double, size_t>> pValuesWithIndex;
    pValuesWithIndex.reserve(clusters.size());

    for (size_t i = 0; i < clusters.size(); ++i) {
        pValuesWithIndex.emplace_back(clusters[i].getPValue(), i);
    }

    std::ranges::sort(pValuesWithIndex);

    std::vector<EvaluatedInteractionCluster> filteredClusters;
    filteredClusters.reserve(clusters.size());

    const auto numPValues = static_cast<double>(pValuesWithIndex.size());
    double minAdjPValue = 1.0;
    for (size_t i = pValuesWithIndex.size(); i-- > 0;) {
        auto rank = static_cast<double>(i + 1);
        double adjustedPValue =
            std::min(minAdjPValue, (pValuesWithIndex[i].first * numPValues / rank));
        minAdjPValue = adjustedPValue;
        clusters[pValuesWithIndex[i].second].setPadj(adjustedPValue);

        if (adjustedPValue <= pAdjThreshold) {
            filteredClusters.emplace_back(std::move(clusters[pValuesWithIndex[i].second]));
        }
    }

    return filteredClusters;
}

}  // namespace pipelines::analyze
