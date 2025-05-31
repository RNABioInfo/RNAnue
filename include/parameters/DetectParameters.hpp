#pragma once

// Standard
#include <cfloat>
#include <cstddef>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "DetectOptions.hpp"
#include "GeneralParameters.hpp"

namespace po = boost::program_options;

namespace pipelines::detect {

struct DetectParameters : public GeneralParameters {
   public:
    size_t maxPrimaryAlignmentCount;
    size_t minimumFragmentLength;
    size_t minimumMapQuality;
    double minimumComplementarity;
    double minimumSiteLengthRatio;
    double maxHybridizationEnergy;
    bool excludeSoftClipping;
    bool removeSplicingEvents;
    bool removeAlternativeSplicing;
    int splicingTolerance;
    bool includeWobbleBasePairsInCrosslinkingSites;
    double minHitGroupContribution;

    DetectParameters(const po::variables_map& params)
        : GeneralParameters(params),
          maxPrimaryAlignmentCount(DetectOptions::maxPrimaryAlignmentCount.extractValue(params)),
          minimumFragmentLength(DetectOptions::minDetectLength.extractValue(params)),
          minimumMapQuality(DetectOptions::minMappingQuality.extractValue(params)),
          minimumComplementarity(DetectOptions::minComplementarity.extractValue(params)),
          minimumSiteLengthRatio(DetectOptions::siteLengthRatio.extractValue(params)),
          maxHybridizationEnergy(DetectOptions::maxEnergy.extractValue(params)),
          excludeSoftClipping(DetectOptions::excludeSoftClipping.extractValue(params)),
          removeSplicingEvents(DetectOptions::filterSplicing.extractValue(params)),
          removeAlternativeSplicing(DetectOptions::allowAltSplicing.extractValue(params)),
          splicingTolerance(DetectOptions::splicingTolerance.extractValue(params)),
          includeWobbleBasePairsInCrosslinkingSites(
              DetectOptions::includeWobble.extractValue(params)),
          minHitGroupContribution(DetectOptions::minHitGroupContribution.extractValue(params)) {};
};

}  // namespace pipelines::detect
