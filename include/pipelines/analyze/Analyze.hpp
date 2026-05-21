#pragma once

// Standard
#include <cstddef>
#include <filesystem>
#include <memory>
#include <vector>

// Boost
#include <boost/math/distributions/binomial.hpp>

// seqan3
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <utility>

// Internal
#include "AnalyzeData.hpp"
#include "AnalyzeParameters.hpp"
#include "AnalyzeSample.hpp"
#include "EvaluatedInteractionCluster.hpp"
#include "FeatureAnnotator.hpp"
#include "TranscriptContributionsByID.hpp"

namespace pipelines::analyze {

using namespace annotation;
namespace math = boost::math;
using seqan3::operator""_tag;

class Analyze {
   public:
    explicit Analyze(AnalyzeParameters params) : parameters(std::move(params)) {};

    void process(const AnalyzeData &data);

   private:
    AnalyzeParameters parameters;
    void processSample(const AnalyzeSample &sample,
                       std::shared_ptr<const FeatureAnnotator> featureAnnotator);

    [[nodiscard]] auto filterByCoverageMetrics(
        std::vector<EvaluatedInteractionCluster> &&clusters) const
        -> std::vector<EvaluatedInteractionCluster>;

    static void parseAnnotatedContiguousFragmentCountsToTranscripts(
        const fs::path &contiguousTranscriptCountsInPath,
        TranscriptContributionsByID &transcriptCounts);

    void parseNonAnnotatedContiguousToSupplementaryFeatures(
        const fs::path &unassignedSingletonsInPath,
        const annotation::FeatureAnnotator &featureAnnotator,
        TranscriptContributionsByID &transcriptCounts);
    static auto parseSampleFragmentCount(const fs::path &sampleCountsInPath) -> float;

    static void writeTranscriptCounts(const TranscriptContributionsByID &featureCounts,
                                      const fs::path &transcriptCountsOutPath);
};

}  // namespace pipelines::analyze
