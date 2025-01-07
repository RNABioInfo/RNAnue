#pragma once
// Standard
#include <cstddef>
#include <string>

namespace constants::pipelines {
const std::string PREPROCESS = "preprocess";
const std::string ALIGN = "align";
const std::string DETECT = "detect";
const std::string ANALYZE = "analyze";
const std::string COMPLETE = "complete";

const std::string GENERAL_DESCRIPTION =
    "RNAnue efficient data analysis for RNA–RNA interactomics.\nRun RNAnue with the subcall "
    "\"complete\" to execute all pipeline steps.\n\nMinimum call: RNAnue complete -t "
    "<treatment-dir> "
    "-o "
    "<output-dir> -f <feature-gff-file> --dbref <reference-genome-file>\nOr run RNAnue with a "
    "config "
    "file: RNAnue complete -c <config-file>\n\nGeneral Options";
const std::string SUBCALL_DESCRIPTION =
    "The subcall to execute. The following subcalls are available: preprocess, align, detect, "
    "analyze, complete.";

const std::string PROCESSING_TREATMENT_MESSAGE = "Processing treatment data";
const std::string PROCESSING_CONTROL_MESSAGE = "Processing control data";

// Preprocess defaults
constexpr size_t defaultChunkSize = 1000000;
constexpr bool trimpolyG = false;
constexpr bool deduplicate = true;
constexpr double defaultAdapterTrimMissmatchRate = 0.05;
constexpr size_t defaultAdapterTrimMinOverlap = 5;
constexpr size_t defaultMinMeanPhreadQuality = 20;
constexpr size_t defaultMinReadLength = 15;
constexpr size_t defaultMinWindowPhredQuality = 20;
constexpr size_t defaultWindowTrimSize = 0;
constexpr size_t defaultMinOverlapMergeSize = 5;
constexpr double defaultMinOverlapMergeMissmatchRate = 0.05;

// Align defaults
constexpr bool defaultMultiMap = false;
constexpr size_t defaultAlignAccuracy = 90;
constexpr size_t defaultMinFragmentScore = 18;
constexpr size_t defaultMinFragmentLength = 20;
constexpr size_t defaultMinSpliceCoverage = 80;

// Detect defaults
constexpr size_t defaultMinMapq = 0;
constexpr double defaultMinComplementarity = 0.5;
constexpr double defaultMinSitelenRatio = 0.1;
constexpr double defaultHybridizationEnergyCutoff = 0;
constexpr int defaultSplicingTolerance = 5;

// Analyze defaults
constexpr double defaultMaxOverlap = 0.5;
constexpr int defaultClusterTolerance = 0;
constexpr double defaultPAdjCutOff = 1.0;
constexpr size_t defaultMinClusterCount = 1;
}  // namespace constants::pipelines

namespace constants::annotation {
constexpr size_t exptectedGffFileTokenCount = 9;
constexpr size_t strandTokenColumn = 6;
}  // namespace constants::annotation
