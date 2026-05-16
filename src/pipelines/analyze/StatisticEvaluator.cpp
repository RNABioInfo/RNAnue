#include "StatisticEvaluator.hpp"

// Standard
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Boost
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/complement.hpp>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "EvaluatedInteractionCluster.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "TranscriptContributionsByID.hpp"
#include "Utility.hpp"

namespace pipelines::analyze {

namespace {

constexpr double backgroundPseudocount = 0.5;

using FeaturePair = std::pair<std::string, std::string>;

struct CandidateContributionSummary {
    double totalContribution{0.0};
    double squaredContributionSum{0.0};

    [[nodiscard]] auto isValid() const noexcept -> bool {
        return std::isfinite(totalContribution) && totalContribution > 0.0 &&
               std::isfinite(squaredContributionSum) && squaredContributionSum > 0.0;
    }

    [[nodiscard]] auto effectiveTrialCount() const noexcept -> double {
        return totalContribution * totalContribution / squaredContributionSum;
    }

    [[nodiscard]] auto effectiveSuccessScale() const noexcept -> double {
        return totalContribution / squaredContributionSum;
    }
};

[[nodiscard]] auto clampProbability(double probability) noexcept -> double {
    if (!std::isfinite(probability)) {
        return 1.0;
    }

    return std::clamp(probability, 0.0, 1.0);
}

[[nodiscard]] auto canonicalFeaturePair(const std::string &firstFeatureID,
                                        const std::string &secondFeatureID) -> FeaturePair {
    if (firstFeatureID <= secondFeatureID) {
        return {firstFeatureID, secondFeatureID};
    }

    return {secondFeatureID, firstFeatureID};
}

[[nodiscard]] auto sanitizedBackgroundContribution(
    const TranscriptContributionsByID &backgroundContributions, const std::string &featureID)
    -> double {
    const auto backgroundIt = backgroundContributions.find(featureID);
    if (backgroundIt == backgroundContributions.end()) {
        Logger::log<LogLevel::DEBUG>(
            "Feature has no contiguous background contribution; using pseudocount. Feature: ",
            featureID);
        return backgroundPseudocount;
    }

    const double contribution = backgroundIt->second;
    if (!std::isfinite(contribution) || contribution < 0.0) {
        Logger::log<LogLevel::WARNING>(
            "Feature background contribution is invalid; using pseudocount. Feature: ", featureID,
            "; contribution: ", contribution);
        return backgroundPseudocount;
    }

    if (helper::isApproxEqual(contribution, 0.0)) {
        Logger::log<LogLevel::DEBUG>(
            "Feature has zero contiguous background contribution; using pseudocount. Feature: ",
            featureID);
    }

    return contribution + backgroundPseudocount;
}

[[nodiscard]] auto summarizeCandidateContributions(
    const std::vector<AnnotatedInteractionCluster> &clusters) -> CandidateContributionSummary {
    CandidateContributionSummary summary{};

    for (const auto &cluster : clusters) {
        const double contribution = cluster.getTranscriptContribution();
        const double squaredContributionSum = cluster.getTranscriptContributionSquaredSum();

        if (!std::isfinite(contribution) || contribution <= 0.0) {
            continue;
        }

        summary.totalContribution += contribution;

        if (std::isfinite(squaredContributionSum) && squaredContributionSum > 0.0) {
            summary.squaredContributionSum += squaredContributionSum;
        } else {
            summary.squaredContributionSum += contribution * contribution;
        }
    }

    return summary;
}

class InteractionNullModel {
   public:
    InteractionNullModel(const std::vector<AnnotatedInteractionCluster> &clusters,
                         const TranscriptContributionsByID &backgroundContributions)
        : clusters(clusters), backgroundContributions(backgroundContributions) {}

    [[nodiscard]] auto clusterProbabilities() const -> std::vector<double> {
        if (clusters.empty()) {
            return {};
        }

        const auto smoothedBackgroundByFeature = buildSmoothedBackgroundByFeature();
        const auto clusterCountsByPair = countClustersByFeaturePair();

        std::vector<double> clusterWeights;
        clusterWeights.reserve(clusters.size());

        double totalWeight = 0.0;
        for (const auto &cluster : clusters) {
            const FeaturePair featurePair =
                canonicalFeaturePair(cluster.getFirstFeatureID(), cluster.getSecondFeatureID());
            const double pairWeight = featurePairWeight(featurePair, smoothedBackgroundByFeature);
            const double clusterCount = static_cast<double>(clusterCountsByPair.at(featurePair));
            const double clusterWeight = pairWeight / clusterCount;

            clusterWeights.emplace_back(clusterWeight);
            totalWeight += clusterWeight;
        }

        if (!std::isfinite(totalWeight) || totalWeight <= 0.0) {
            Logger::log<LogLevel::WARNING>(
                "Interaction null model has no valid background weight; using uniform cluster "
                "probabilities.");
            return std::vector<double>(clusters.size(), 1.0 / static_cast<double>(clusters.size()));
        }

        std::vector<double> probabilities;
        probabilities.reserve(clusterWeights.size());
        for (const double clusterWeight : clusterWeights) {
            probabilities.emplace_back(clampProbability(clusterWeight / totalWeight));
        }

        return probabilities;
    }

   private:
    const std::vector<AnnotatedInteractionCluster> &clusters;
    const TranscriptContributionsByID &backgroundContributions;

    [[nodiscard]] auto buildSmoothedBackgroundByFeature() const
        -> std::unordered_map<std::string, double> {
        std::unordered_map<std::string, double> smoothedBackgroundByFeature;

        for (const auto &cluster : clusters) {
            addFeatureBackground(cluster.getFirstFeatureID(), smoothedBackgroundByFeature);
            addFeatureBackground(cluster.getSecondFeatureID(), smoothedBackgroundByFeature);
        }

        return smoothedBackgroundByFeature;
    }

    void addFeatureBackground(
        const std::string &featureID,
        std::unordered_map<std::string, double> &smoothedBackgroundByFeature) const {
        if (smoothedBackgroundByFeature.contains(featureID)) {
            return;
        }

        smoothedBackgroundByFeature.emplace(
            featureID, sanitizedBackgroundContribution(backgroundContributions, featureID));
    }

    [[nodiscard]] auto countClustersByFeaturePair() const -> std::map<FeaturePair, size_t> {
        std::map<FeaturePair, size_t> clusterCountsByPair;

        for (const auto &cluster : clusters) {
            ++clusterCountsByPair[canonicalFeaturePair(cluster.getFirstFeatureID(),
                                                       cluster.getSecondFeatureID())];
        }

        return clusterCountsByPair;
    }

    [[nodiscard]] static auto featurePairWeight(
        const FeaturePair &featurePair,
        const std::unordered_map<std::string, double> &smoothedBackgroundByFeature) -> double {
        const double firstFeatureWeight = smoothedBackgroundByFeature.at(featurePair.first);
        const double secondFeatureWeight = smoothedBackgroundByFeature.at(featurePair.second);

        if (featurePair.first == featurePair.second) {
            return firstFeatureWeight * secondFeatureWeight;
        }

        return 2.0 * firstFeatureWeight * secondFeatureWeight;
    }
};

[[nodiscard]] auto binomialUpperTail(double trials, double probability,
                                     double observedEffectiveSuccesses) -> double {
    if (!std::isfinite(trials) || trials <= 0.0 || !std::isfinite(probability) ||
        !std::isfinite(observedEffectiveSuccesses)) {
        return 1.0;
    }

    probability = std::clamp(probability, 0.0, 1.0);

    if (observedEffectiveSuccesses <= 0.0) {
        return 1.0;
    }

    if (helper::isApproxEqual(probability, 0.0)) {
        return 0.0;
    }

    if (helper::isApproxEqual(probability, 1.0)) {
        return observedEffectiveSuccesses <= trials ? 1.0 : 0.0;
    }

    const auto distribution = math::binomial_distribution<double>(trials, probability);
    const double threshold = observedEffectiveSuccesses - 1.0;

    if (threshold < 0.0) {
        return 1.0;
    }

    if (threshold >= trials) {
        return 0.0;
    }

    return clampProbability(math::cdf(math::complement(distribution, threshold)));
}

}  // namespace

auto StatisticEvaluator::evaluate(std::vector<AnnotatedInteractionCluster> &clusters,
                                  const TranscriptContributionsByID &backgroundContributions,
                                  double padjThreshold)
    -> std::vector<EvaluatedInteractionCluster> {
    Logger::log("Evaluating ", clusters.size(),
                " clusters with abundance-corrected p-adj threshold: ", padjThreshold);
    auto evaluatedClusters = evaluatePValues(clusters, backgroundContributions);
    return evaluatePAdjValues(evaluatedClusters, padjThreshold);
}

auto StatisticEvaluator::evaluatePValues(std::vector<AnnotatedInteractionCluster> &clusters,
                                         const TranscriptContributionsByID &backgroundContributions)
    -> std::vector<EvaluatedInteractionCluster> {
    std::vector<EvaluatedInteractionCluster> evaluatedClusters;
    evaluatedClusters.reserve(clusters.size());

    if (clusters.empty()) {
        return evaluatedClusters;
    }

    const CandidateContributionSummary contributionSummary =
        summarizeCandidateContributions(clusters);
    const std::vector<double> clusterProbabilities =
        InteractionNullModel{clusters, backgroundContributions}.clusterProbabilities();

    if (!contributionSummary.isValid()) {
        Logger::log<LogLevel::WARNING>(
            "No valid candidate interaction contribution available for p-value calculation; "
            "assigning p=1 to all clusters.");
        for (auto &cluster : clusters) {
            evaluatedClusters.emplace_back(std::move(cluster), 1.0);
        }
        return evaluatedClusters;
    }

    const double effectiveTrialCount = contributionSummary.effectiveTrialCount();
    const double effectiveSuccessScale = contributionSummary.effectiveSuccessScale();

    for (size_t i = 0; i < clusters.size(); ++i) {
        auto &cluster = clusters[i];
        const double observedContribution = cluster.getTranscriptContribution();
        const double pValue =
            (!std::isfinite(observedContribution) || observedContribution <= 0.0)
                ? 1.0
                : binomialUpperTail(effectiveTrialCount, clusterProbabilities[i],
                                    observedContribution * effectiveSuccessScale);
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
        adjustedPValue = clampProbability(adjustedPValue);
        minAdjPValue = adjustedPValue;
        clusters[pValuesWithIndex[i].second].setPadj(adjustedPValue);
    }

    for (auto &cluster : clusters) {
        if (cluster.getPadj() <= pAdjThreshold) {
            filteredClusters.emplace_back(std::move(cluster));
        }
    }

    return filteredClusters;
}

}  // namespace pipelines::analyze
