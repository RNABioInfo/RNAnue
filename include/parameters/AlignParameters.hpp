#pragma once

// Standard
#include <cstddef>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "AlignOptions.hpp"
#include "GeneralParameters.hpp"

namespace po = boost::program_options;

namespace pipelines::align {

enum class AlignmentBackend { Segemehl, Star };

struct AlignParameters : public GeneralParameters {
    std::filesystem::path referenceGenome;
    bool multimapAlignments;
    size_t minLengthThreshold;
    size_t accuracy;
    size_t minimumFragmentScore;
    size_t minimumFragmentLength;
    size_t minimumSpliceCoverage;
    AlignmentBackend aligner;
    int starMaxMultimaps;
    std::optional<int> starMinJunctionOverhang;
    int starMaxSegmentGap;
    int starMinNonchimericScoreDrop;
    int starMaxChimericScoreDrop;
    int starMaxIntronLength;

    AlignParameters(const po::variables_map& params)
        : GeneralParameters(params),
          referenceGenome(AlignOptions::refGenome.extractValue(params)),
          multimapAlignments(AlignOptions::allowMultimap.extractValue(params)),
          minLengthThreshold(AlignOptions::minAlignLength.extractValue(params)),
          accuracy(AlignOptions::accuracy.extractValue(params)),
          minimumFragmentScore(AlignOptions::minFragmentScore.extractValue(params)),
          minimumFragmentLength(AlignOptions::minFragmentLength.extractValue(params)),
          minimumSpliceCoverage(AlignOptions::minSpliceCoverage.extractValue(params)),
          aligner(validateAligner(AlignOptions::aligner.extractValue(params))),
          starMaxMultimaps(AlignOptions::starMaxMultimaps.extractValue(params)),
          starMinJunctionOverhang(validateJunctionOverhang(
              AlignOptions::starMinJunctionOverhang.extractValue(params))),
          starMaxSegmentGap(AlignOptions::starMaxSegmentGap.extractValue(params)),
          starMinNonchimericScoreDrop(AlignOptions::starMinNonchimericScoreDrop.extractValue(params)),
          starMaxChimericScoreDrop(AlignOptions::starMaxChimericScoreDrop.extractValue(params)),
          starMaxIntronLength(AlignOptions::starMaxIntronLength.extractValue(params)) {
        if (aligner == AlignmentBackend::Star &&
            (minimumFragmentLength == 0 || minimumFragmentLength > INT_MAX)) {
            throw std::invalid_argument("min_fragment_length must be between 1 and INT_MAX for STAR");
        }
        std::apply([&](const auto&... option) {
            (recordExplicitStarOption(params, option.getLongName()), ...);
        }, AlignOptions::starOptions);
    }

    // Keep configured values intact; these are the intended future STAR values.
    [[nodiscard]] auto effectiveStarMaxMultimaps() const -> int {
        return multimapAlignments ? starMaxMultimaps : 1;
    }

    [[nodiscard]] auto effectiveStarMinJunctionOverhang() const -> size_t {
        return starMinJunctionOverhang ? static_cast<size_t>(*starMinJunctionOverhang)
                                      : minimumFragmentLength;
    }

    // Called only by alignment-bearing pipelines, before any pipeline side effects.
    // Extraction itself remains usable without an operational STAR backend.
    void validateBackendAvailability() const {
        if (aligner == AlignmentBackend::Star) {
            throw std::invalid_argument("STAR alignment backend is not implemented yet; use --aligner segemehl");
        }
        if (explicitStarOption) {
            throw std::invalid_argument("--" + *explicitStarOption +
                                        " requires aligner=star; STAR controls cannot be used with aligner=segemehl");
        }
    }

   private:
    std::optional<std::string> explicitStarOption;

    static auto validateAligner(const std::string& value) -> AlignmentBackend {
        if (value == "segemehl") return AlignmentBackend::Segemehl;
        if (value == "star") return AlignmentBackend::Star;
        throw std::invalid_argument("aligner must be exactly 'segemehl' or 'star'");
    }

    static auto validateJunctionOverhang(std::optional<int> value) -> std::optional<int> {
        if (value && *value < 1) {
            throw std::invalid_argument("star_min_junction_overhang must be between 1 and INT_MAX");
        }
        return value;
    }

    void recordExplicitStarOption(const po::variables_map& params, const std::string& name) {
        const auto entry = params.find(name);
        if (!explicitStarOption && entry != params.end() && !entry->second.defaulted()) {
            explicitStarOption = name;
        }
    }
};

}  // namespace pipelines::align
