#pragma once

// Standard
#include <cstddef>
#include <filesystem>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "AlignOptions.hpp"
#include "GeneralParameters.hpp"

namespace po = boost::program_options;

namespace pipelines::align {

struct AlignParameters : public GeneralParameters {
    std::filesystem::path referenceGenome;
    bool multimapAlignments;
    size_t minLengthThreshold;
    size_t accuracy;
    size_t minimumFragmentScore;
    size_t minimumFragmentLength;
    size_t minimumSpliceCoverage;

    AlignParameters(const po::variables_map& params)
        : GeneralParameters(params),
          referenceGenome(AlignOptions::refGenome.extractValue(params)),
          multimapAlignments(AlignOptions::allowMultimap.extractValue(params)),
          minLengthThreshold(AlignOptions::minAlignLength.extractValue(params)),
          accuracy(AlignOptions::accuracy.extractValue(params)),
          minimumFragmentScore(AlignOptions::minFragmentScore.extractValue(params)),
          minimumFragmentLength(AlignOptions::minFragmentLength.extractValue(params)),
          minimumSpliceCoverage(AlignOptions::minSpliceCoverage.extractValue(params)) {};
};

}  // namespace pipelines::align
