#include "StatisticEvaluator.hpp"

#include "Logger.hpp"

namespace pipelines::analyze {

auto StatisticEvaluator::evaluate(std::vector<AnnotatedInteractionCluster> &clusters,
                                  const std::unordered_map<std::string, size_t> &transcriptCounts,
                                  size_t totalTranscriptCount, double padjThreshold)
    -> std::vector<EvaluatedInteractionCluster> {
    auto evaluatedClusters = evaluatePValues(clusters, transcriptCounts, totalTranscriptCount);
    return evaluatePAdjValues(evaluatedClusters, padjThreshold);
}

auto StatisticEvaluator::getTranscriptProbabilities(
    const std::unordered_map<std::string, size_t> &transcriptCounts,
    const size_t totalTranscriptCount) -> std::unordered_map<std::string, double> {
    std::unordered_map<std::string, double> transcriptProbabilities;

    double cumulativeProbability = 0;
    for (const auto &[transcriptID, count] : transcriptCounts) {
        transcriptProbabilities[transcriptID] =
            static_cast<double>(count) / static_cast<double>(totalTranscriptCount);
        cumulativeProbability += transcriptProbabilities[transcriptID];
    }

    // Normalize probabilities to sum to 1
    for (auto &[transcriptID, probability] : transcriptProbabilities) {
        probability /= cumulativeProbability;
    }

    return transcriptProbabilities;
}

auto StatisticEvaluator::evaluatePValues(
    std::vector<AnnotatedInteractionCluster> &clusters,
    const std::unordered_map<std::string, size_t> &transcriptCounts,
    const size_t totalTranscriptCount) -> std::vector<EvaluatedInteractionCluster> {
    const auto transcriptProbabilities =
        getTranscriptProbabilities(transcriptCounts, totalTranscriptCount);

    std::vector<EvaluatedInteractionCluster> evaluatedClusters;
    evaluatedClusters.reserve(clusters.size());

    for (auto &cluster : clusters) {
        const auto &firstTranscriptID = cluster.getFirstFeatureID();
        const auto &secondTranscriptID = cluster.getSecondFeatureID();

        auto findProbability = [&](const std::string &transcriptID) -> std::optional<double> {
            auto iterator = transcriptProbabilities.find(transcriptID);
            return (iterator != transcriptProbabilities.end())
                       ? std::optional<double>{iterator->second}
                       : std::nullopt;
        };

        auto firstTranscriptProbability = findProbability(firstTranscriptID);
        auto secondTranscriptProbability = findProbability(secondTranscriptID);

        if (!firstTranscriptProbability || !secondTranscriptProbability) {
            Logger::log(LogLevel::WARNING,
                        "Could not find transcript probabilities for cluster with transcripts: ",
                        firstTranscriptID, ", ", secondTranscriptID);
            continue;
        }

        const double ligationByChanceProbability =
            (firstTranscriptID == secondTranscriptID)
                ? firstTranscriptProbability.value() * secondTranscriptProbability.value()
                : 2 * firstTranscriptProbability.value() * secondTranscriptProbability.value();
        ;

        const auto binomialDistribution =
            math::binomial_distribution((double)totalTranscriptCount, ligationByChanceProbability);

        const double pValue = 1 - math::cdf(binomialDistribution, cluster.fragmentCount());

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
