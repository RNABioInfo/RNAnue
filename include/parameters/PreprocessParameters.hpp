#pragma once

// Standard
#include <cstddef>
#include <filesystem>
#include <string>
#include <variant>
#include <vector>

// Boost
#include <boost/program_options/variables_map.hpp>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/utility/range/to.hpp>

// Internal
#include "GeneralParameters.hpp"
#include "PreprocessOptions.hpp"

namespace pipelines::preprocess {

namespace po = boost::program_options;

using AdapterInput = std::variant<std::monostate, std::filesystem::path, seqan3::dna5_vector>;

class PreprocessParameters : public GeneralParameters {
   public:
    bool preprocessEnabled;

    bool trimPolyG;
    bool deduplicate;

    AdapterInput adapter5Forward;
    AdapterInput adapter3Forward;
    AdapterInput adapter5Reverse;
    AdapterInput adapter3Reverse;

    size_t minPolyGCount;
    double maxMissMatchFractionTrimming;
    size_t minOverlapTrimming;
    size_t minQualityThreshold;
    size_t minLengthThreshold;
    size_t minMeanWindowQuality;
    size_t windowTrimmingSize;
    size_t minOverlapMerging;
    double maxMissMatchFractionMerging;

    PreprocessParameters(const po::variables_map& params)
        : GeneralParameters(params),
          preprocessEnabled(PreprocessOptions::enablePreprocess.extractValue(params)),
          trimPolyG(PreprocessOptions::trimPolyG.extractValue(params)),
          deduplicate(PreprocessOptions::enableDeduplicate.extractValue(params)),
          adapter5Forward(getAdapter(PreprocessOptions::adpt5f.extractValue(params))),
          adapter3Forward(getAdapter(PreprocessOptions::adpt3f.extractValue(params))),
          adapter5Reverse(getAdapter(PreprocessOptions::adpt5r.extractValue(params))),
          adapter3Reverse(getAdapter(PreprocessOptions::adpt3r.extractValue(params))),
          minPolyGCount(PreprocessOptions::minPolyGCount.extractValue(params)),
          maxMissMatchFractionTrimming(PreprocessOptions::mtrim.extractValue(params)),
          minOverlapTrimming(PreprocessOptions::minOvlTrim.extractValue(params)),
          minQualityThreshold(PreprocessOptions::minQual.extractValue(params)),
          minLengthThreshold(PreprocessOptions::minLen.extractValue(params)),
          minMeanWindowQuality(PreprocessOptions::wqual.extractValue(params)),
          windowTrimmingSize(PreprocessOptions::wtrim.extractValue(params)),
          minOverlapMerging(PreprocessOptions::minOvl.extractValue(params)),
          maxMissMatchFractionMerging(PreprocessOptions::mmerge.extractValue(params)) {}

   private:
    static auto getAdapter(const std::string& adapterStr) -> AdapterInput {
        if (std::filesystem::exists(adapterStr)) {
            return std::filesystem::path(adapterStr);
        }

        return adapterStr | seqan3::views::char_to<seqan3::dna5> |
               seqan3::ranges::to<std::vector>();
    }
};

}  // namespace pipelines::preprocess
