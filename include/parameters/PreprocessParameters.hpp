#pragma once

// Standard
#include <cstddef>
#include <cstdint>
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
#include "Constants.hpp"
#include "GeneralParameters.hpp"
#include "ParameterValidator.hpp"

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
          preprocessEnabled(validatePreprocessEnabled(params)),
          trimPolyG(validateTrimPolyG(params)),
          deduplicate(validateDedupliate(params)),
          adapter5Forward(validateAdapter(params, "adpt5f")),
          adapter3Forward(validateAdapter(params, "adpt3f")),
          adapter5Reverse(validateAdapter(params, "adpt5r")),
          adapter3Reverse(validateAdapter(params, "adpt3r")),
          minPolyGCount(
              ParameterValidator::validateArithmetic(params, "minpolygcount", size_t{0}, SIZE_MAX)),
          maxMissMatchFractionTrimming(
              ParameterValidator::validateArithmetic(params, "mtrim", 0.0, 1.0)),
          minOverlapTrimming(
              ParameterValidator::validateArithmetic(params, "minovltrim", size_t{0}, SIZE_MAX)),
          minQualityThreshold(
              ParameterValidator::validateArithmetic(params, "minqual", size_t{0}, SIZE_MAX)),
          minLengthThreshold(ParameterValidator::validateArithmetic<size_t>(params, "minlen",
                                                                            size_t{1}, SIZE_MAX)),
          minMeanWindowQuality(
              ParameterValidator::validateArithmetic(params, "wqual", size_t{0}, SIZE_MAX)),
          windowTrimmingSize(
              ParameterValidator::validateArithmetic(params, "wtrim", size_t{0}, SIZE_MAX)),
          minOverlapMerging(
              ParameterValidator::validateArithmetic(params, "minovl", size_t{0}, SIZE_MAX)),
          maxMissMatchFractionMerging(
              ParameterValidator::validateArithmetic(params, "mmerge", 0.0, 1.0)) {}

   private:
    static auto validateTrimPolyG(const po::variables_map& params) -> bool {
        return params["trimpolyg"].as<bool>();
    }

    static auto validateDedupliate(const po::variables_map& params) -> bool {
        return params["deduplicate"].as<bool>();
    }

    static auto validatePreprocessEnabled(const po::variables_map& params) -> bool {
        return params[constants::pipelines::PREPROCESS].as<bool>();
    }

    static auto validateAdapter(const po::variables_map& params, const std::string& paramName)
        -> AdapterInput {
        if (params.count(paramName) != 0U) {
            const std::string adapterStr = params[paramName].as<std::string>();

            if (std::filesystem::exists(adapterStr)) {
                return std::filesystem::path(adapterStr);
            }

            return adapterStr | seqan3::views::char_to<seqan3::dna5> |
                   seqan3::ranges::to<std::vector>();
        }
        return std::monostate{};
    }
};

}  // namespace pipelines::preprocess
