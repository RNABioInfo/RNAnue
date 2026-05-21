#pragma once
// Standard
#include <cstddef>
#include <string>
#include <string_view>
#include <unordered_set>

namespace constants::pipelines {
const std::string SUBCALL_PARAMETER_KEY = "subcall";

const std::string PREPROCESS = "preprocess";
const std::string ALIGN = "align";
const std::string DETECT = "detect";
const std::string ANALYZE = "analyze";
const std::string POSTPROCESS = "postprocess";
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
constexpr bool defaultTrimpolyG = false;
constexpr size_t dfaultPolyGCutoff = 5;
constexpr bool defaultDeduplicate = true;
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
constexpr std::string_view maskedAnnotationFileName = "masked_annotations.gtf";
constexpr std::string_view maskedReferenceGenomeFileName = "masked_reference_genome.fa";

// Detect defaults
constexpr size_t defaultMinMapq = 0;
constexpr double defaultMinComplementarity = 0.5;
constexpr double defaultMinSitelenRatio = 0.1;
constexpr double defaultHybridizationEnergyCutoff = 0;
constexpr int defaultSplicingTolerance = 5;

// Analyze defaults
constexpr float clusterOverlapFractionMin = 1e-9;
constexpr double defaultMaxOverlap = 0.1;
constexpr int defaultClusterTolerance = 0;
constexpr double defaultPAdjCutOff = 1.0;
constexpr size_t defaultMinClusterCount = 1;
}  // namespace constants::pipelines

namespace constants::annotation {
inline static const std::unordered_set<std::string> defaultAllTranscriptTypes{"antisense_RNA",
                                                                              "C_gene_segment",
                                                                              "D_gene_segment",
                                                                              "exon",
                                                                              "five_prime_UTR",
                                                                              "telomerase_RNA",
                                                                              "three_prime_UTR",
                                                                              "gene",
                                                                              "J_gene_segment",
                                                                              "lnc_RNA",
                                                                              "miRNA",
                                                                              "mRNA",
                                                                              "primary_transcript",
                                                                              "pseudogene",
                                                                              "RNase_MRP_RNA",
                                                                              "RNase_P_RNA",
                                                                              "rRNA",
                                                                              "scaRNA",
                                                                              "scRNA",
                                                                              "snoRNA",
                                                                              "snRNA",
                                                                              "transcript",
                                                                              "tmRNA",
                                                                              "tRNA",
                                                                              "vault_RNA",
                                                                              "V_gene_segment",
                                                                              "Y_RNA",
                                                                              "TU",
                                                                              "transcription_unit"};
constexpr size_t expectedAnnotationFileTokenCount = 9;
constexpr size_t strandTokenColumn = 6;
}  // namespace constants::annotation

namespace constants::interaction {
static const std::string clusterIDHeader = "cluster_ID";
static const std::string firstFeatureIDHeader = "fst_feat_id";
static const std::string firstSegmentReferenceHeader = "fst_seg_chr";
static const std::string firstSegmentStartHeader = "fst_seg_strt";
static const std::string firstSegmentEndHeader = "fst_seg_end";
static const std::string firstSegmentStrandHeader = "fst_seg_strd";
static const std::string secondFeatureIDHeader = "sec_feat_id";
static const std::string secondSegmentReferenceHeader = "sec_seg_chr";
static const std::string secondSegmentStartHeader = "sec_seg_strt";
static const std::string secondSegmentEndHeader = "sec_seg_end";
static const std::string secondSegmentStrandHeader = "sec_seg_strd";
static const std::string transcriptContributionHeader = "no_splits";
static const std::string totalSpanBpHeader = "total_span_bp";
static const std::string effectiveCoverageSpanBpHeader = "effective_coverage_span_bp";
static const std::string supportPerTotalBpHeader = "support_per_total_bp";
static const std::string supportPerEffectiveBpHeader = "support_per_effective_bp";
static const std::string coverageConcentrationHeader = "coverage_concentration";
static const std::string coverageComponentsHeader = "coverage_components";
static const std::string armBalanceHeader = "arm_balance";
static const std::string coverageProfileHeader = "coverage_profile";
static const std::string meanInterCrosslinkCountHeader = "mean_inter_crosslinks";
static const std::string sdInterCrosslinksHeader = "sd_inter_crosslinks";
static const std::string globalComplementarityScoreHeader = "gcs";
static const std::string globalHybridizationScoreHeader = "ghs";
static const std::string pValueHeader = "p_value";
static const std::string padjValueHeader = "padj_value";
}  // namespace constants::interaction
