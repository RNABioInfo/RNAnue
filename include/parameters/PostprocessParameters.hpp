#pragma once

// Internal
#include <boost/program_options/variables_map.hpp>

#include "GeneralParameters.hpp"
#include "PostprocessOptions.hpp"

namespace pipelines::postprocess {

namespace po = boost::program_options;

class PostprocessParameters : public GeneralParameters {
   public:
    PostprocessParameters(const po::variables_map& params)
        : GeneralParameters(params),
          minSegmentOverlapFraction(
              PostprocessOptions::minSegmentFractionOverlap.extractValue(params)) {}

    float minSegmentOverlapFraction;
};

}  // namespace pipelines::postprocess
