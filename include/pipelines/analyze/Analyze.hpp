#pragma once

// Standard
#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>

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
#include "FeatureAnnotator.hpp"

namespace pipelines::analyze {

namespace math = boost::math;
using seqan3::operator""_tag;
using namespace annotation;

class Analyze {
   public:
    explicit Analyze(AnalyzeParameters params) : parameters(std::move(params)) {};

    void process(const AnalyzeData &data);

   private:
    AnalyzeParameters parameters;

    void processSample(const AnalyzeSample &sample,
                       std::shared_ptr<const FeatureAnnotator> featureAnnotator);

    static void assignAnnotatedContiguousFragmentCountsToTranscripts(
        const fs::path &contiguousTranscriptCountsInPath,
        std::unordered_map<std::string, size_t> &transcriptCounts);

    void assignNonAnnotatedContiguousToSupplementaryFeatures(
        const fs::path &unassignedSingletonsInPath, annotation::FeatureAnnotator &featureAnnotator,
        std::unordered_map<std::string, size_t> &transcriptCounts);
    static auto parseSampleFragmentCount(const fs::path &sampleCountsInPath) -> size_t;

    static void writeTranscriptCounts(const std::unordered_map<std::string, size_t> &featureCounts,
                                      const fs::path &transcriptCountsOutPath);
};

}  // namespace pipelines::analyze
