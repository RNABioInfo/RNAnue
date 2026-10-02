#pragma once

// Standard
#include <unistd.h>

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
#include "AlignmentBackend.hpp"
#include "GeneralParameters.hpp"

namespace po = boost::program_options;

namespace pipelines::align {

struct AlignParameters : public GeneralParameters {
    std::filesystem::path referenceGenome;
    std::optional<std::filesystem::path> alignmentIndex;
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
          alignmentIndex(AlignOptions::alignmentIndex.extractValue(params)),
          multimapAlignments(AlignOptions::allowMultimap.extractValue(params)),
          minLengthThreshold(AlignOptions::minAlignLength.extractValue(params)),
          accuracy(AlignOptions::accuracy.extractValue(params)),
          minimumFragmentScore(AlignOptions::minFragmentScore.extractValue(params)),
          minimumFragmentLength(AlignOptions::minFragmentLength.extractValue(params)),
          minimumSpliceCoverage(AlignOptions::minSpliceCoverage.extractValue(params)),
          aligner(validateAligner(AlignOptions::aligner.extractValue(params))),
          starMaxMultimaps(AlignOptions::starMaxMultimaps.extractValue(params)),
          starMinJunctionOverhang(
              validateJunctionOverhang(AlignOptions::starMinJunctionOverhang.extractValue(params))),
          starMaxSegmentGap(AlignOptions::starMaxSegmentGap.extractValue(params)),
          starMinNonchimericScoreDrop(
              AlignOptions::starMinNonchimericScoreDrop.extractValue(params)),
          starMaxChimericScoreDrop(AlignOptions::starMaxChimericScoreDrop.extractValue(params)),
          starMaxIntronLength(AlignOptions::starMaxIntronLength.extractValue(params)) {
        if (alignmentIndex && maskMultiCopyGenes) {
            throw std::invalid_argument(
                "--alignment_index is incompatible with mask_multicopy_genes=true; "
                "use --mask_multicopy_genes=false or mask_multicopy_genes = false in configuration");
        }

        if (alignmentIndex && alignmentIndex->empty()) {
            throw std::invalid_argument("--alignment_index must not be empty");
        }

        if (aligner == AlignmentBackend::Star &&
            (minimumFragmentLength == 0 || minimumFragmentLength > INT_MAX)) {
            throw std::invalid_argument(
                "min_fragment_length must be between 1 and INT_MAX for STAR");
        }
        std::apply(
            [&](const auto&... option) {
                (recordExplicitStarOption(params, option.getLongName()), ...);
            },
            AlignOptions::starOptions);
        for (const auto& name :
             {AlignOptions::accuracy.getLongName(), AlignOptions::minFragmentScore.getLongName()}) {
            const auto entry = params.find(name);
            if (entry != params.end() && !entry->second.defaulted()) explicitSegemehlOption = name;
        }
    }

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
        if (alignmentIndex) {
            std::error_code error;
            const bool isStar = aligner == AlignmentBackend::Star;
            const bool correctType = isStar
                                         ? std::filesystem::is_directory(*alignmentIndex, error)
                                         : std::filesystem::is_regular_file(*alignmentIndex, error);
            if (!correctType || error ||
                access(alignmentIndex->c_str(), isStar ? R_OK | X_OK : R_OK) != 0) {
                throw std::invalid_argument(
                    "--alignment_index must be an accessible " +
                    std::string(isStar ? "STAR index directory: " : "segemehl index file: ") +
                    alignmentIndex->string());
            }
        }
        if (aligner == AlignmentBackend::Star) {
            if (explicitSegemehlOption) {
                throw std::invalid_argument("--" + *explicitSegemehlOption +
                                            " requires aligner=segemehl");
            }
            return;
        }
        if (explicitStarOption) {
            throw std::invalid_argument(
                "--" + *explicitStarOption +
                " requires aligner=star; STAR controls cannot be used with aligner=segemehl");
        }
    }

   private:
    std::optional<std::string> explicitStarOption;
    std::optional<std::string> explicitSegemehlOption;

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
