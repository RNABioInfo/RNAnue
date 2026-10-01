#include <gtest/gtest.h>

#include <climits>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "CompleteParameters.hpp"
#include "OtherOptions.hpp"
#include "ParameterOptions.hpp"
#include "SubcallOptions.hpp"
#include "Utility.hpp"

namespace {
using pipelines::align::AlignParameters;
using pipelines::align::AlignmentBackend;

const auto parameterOptions = std::tuple_cat(
    GeneralOptions::allOptions, PreprocessOptions::allOptions, AlignOptions::allOptions,
    DetectOptions::allOptions, AnalyzeOptions::allOptions, PostprocessOptions::allOptions);

auto options() -> po::options_description {
    po::options_description result;
    std::apply([&](const auto&... option) { (option.addOptionTo(result), ...); }, parameterOptions);
    result.add(ParameterOptions::getOptions<OtherOptions>())
          .add(ParameterOptions::getOptions<SubcallOptions>());
    return result;
}

class ParameterOptionsTests : public testing::Test {
   protected:
    std::filesystem::path root = std::filesystem::temp_directory_path() /
                                ("RNAnue-parameter-tests-" + helper::getUUID());
    void SetUp() override {
        std::filesystem::create_directories(root);
        std::ofstream(root / "features.gff") << "##gff-version 3\n";
    }
    void TearDown() override { std::filesystem::remove_all(root); }

    auto parse(std::vector<std::string> args = {}, const std::string& config = "")
        -> po::variables_map {
        // Required general parameters are valid paths; output need not exist.
        const std::string required = "treatment_dir = " + root.string() +
            "\noutput_dir = " + (root / "out").string() +
            "\nfeatures = " + (root / "features.gff").string() +
            "\nreference_genome = reference.fa\n";
        auto description = options();
        po::variables_map result;
        po::store(po::command_line_parser(args).options(description)
                      .style(po::command_line_style::default_style & ~po::command_line_style::allow_guessing)
                      .run(), result);
        po::notify(result);
        std::istringstream input(config);
        po::store(po::parse_config_file(input, description), result);
        std::istringstream fallback(required);
        po::store(po::parse_config_file(fallback, description), result);
        po::notify(result);
        return result;
    }
};

TEST_F(ParameterOptionsTests, CanonicalNamesAndShortcutsAreUnique) {
    auto registered = options();
    std::set<std::string> names;
    for (const auto& option : registered.options()) {
        EXPECT_TRUE(std::regex_match(option->long_name(), std::regex("[a-z][a-z0-9]*(_[a-z0-9]+)*")))
            << option->long_name();
        EXPECT_TRUE(names.insert(option->long_name()).second) << option->long_name();
    }
    std::map<char, std::string> shorts;
    const auto all = std::tuple_cat(parameterOptions, OtherOptions::allOptions, SubcallOptions::allOptions);
    std::apply([&](const auto&... option) {
        const auto check = [&](const auto& value) {
            if (const auto shortName = value.getShortName()) {
                EXPECT_TRUE(shorts.emplace(*shortName, value.getLongName()).second);
            }
        };
        (check(option), ...);
    }, all);
    const std::map<char, std::string> expected{
        {'T', "treatment_dir"}, {'C', "control_dir"}, {'t', "threads"}, {'o', "output_dir"},
        {'r', "reference_genome"}, {'f', "features"}, {'a', "aligner"}, {'i', "alignment_index"}, {'c', "config"},
        {'q', "min_read_quality"}, {'l', "min_read_length"}, {'h', "help"}, {'v', "version"}};
    EXPECT_EQ(shorts, expected);
}

TEST_F(ParameterOptionsTests, EveryOptionAppearsInHelpWithItsDescriptionAndDefault) {
    std::ostringstream help;
    help << options();
    std::apply([&](const auto&... option) {
        const auto check = [&](const auto& value) {
            EXPECT_NE(help.str().find("--" + value.getLongName()), std::string::npos);
            if constexpr (requires { value.getDefaultValue(); }) {
                EXPECT_NE(value.getDescription().find("default: "), std::string::npos);
            }
        };
        (check(option), ...);
    }, parameterOptions);
}

TEST_F(ParameterOptionsTests, StarDefaultsDoNotBlockSegemehl) {
    AlignParameters p{parse()};
    EXPECT_EQ(p.aligner, AlignmentBackend::Segemehl);
    EXPECT_EQ(p.starMaxMultimaps, 10);
    EXPECT_FALSE(p.starMinJunctionOverhang);
    EXPECT_EQ(p.effectiveStarMinJunctionOverhang(), 15);
    EXPECT_EQ(p.starMaxSegmentGap, 3);
    EXPECT_EQ(p.starMinNonchimericScoreDrop, 10);
    EXPECT_EQ(p.starMaxChimericScoreDrop, 30);
    EXPECT_EQ(p.starMaxIntronLength, 10);
    EXPECT_NO_THROW(p.validateBackendAvailability());
}

TEST_F(ParameterOptionsTests, StarCanBeExtractedIndependentlyOfExecution) {
    AlignParameters p{parse({"--aligner=star", "--star_max_multimaps=17", "--min_fragment_length=23",
                            "--star_max_segment_gap=4", "--star_min_nonchimeric_score_drop=8",
                            "--star_max_chimeric_score_drop=26", "--star_max_intron_length=1000"})};
    EXPECT_EQ(p.aligner, AlignmentBackend::Star);
    EXPECT_EQ(p.starMaxMultimaps, 17);
    EXPECT_EQ(p.effectiveStarMaxMultimaps(), 17);
    EXPECT_EQ(p.effectiveStarMinJunctionOverhang(), 23);
    EXPECT_EQ(p.starMaxSegmentGap, 4);
    EXPECT_EQ(p.starMinNonchimericScoreDrop, 8);
    EXPECT_EQ(p.starMaxChimericScoreDrop, 26);
    EXPECT_EQ(p.starMaxIntronLength, 1000);
    EXPECT_NO_THROW(p.validateBackendAvailability());
}

TEST_F(ParameterOptionsTests, OverhangOverrideAndMultimapCapRemainConfigured) {
    AlignParameters p{parse({"--aligner=star", "--min_fragment_length=23",
                            "--star_min_junction_overhang=9", "--star_max_multimaps=17",
                            "--no_multimapping"})};
    EXPECT_EQ(p.starMinJunctionOverhang, 9);
    EXPECT_EQ(p.effectiveStarMinJunctionOverhang(), 9);
    EXPECT_EQ(p.starMaxMultimaps, 17);
    EXPECT_EQ(p.effectiveStarMaxMultimaps(), 1);
}

TEST_F(ParameterOptionsTests, StarIntegerBoundaries) {
    std::apply([&](const auto&... option) {
        const auto check = [&](const auto& value) {
            const auto name = value.getLongName();
            SCOPED_TRACE(name);
            AlignParameters p{parse({"--aligner=star", "--" + name + "=" + std::to_string(INT_MAX)})};
            EXPECT_EQ(value.extractValue(parse({"--" + name + "=" + std::to_string(INT_MAX)})), INT_MAX);
            EXPECT_THROW(parse({"--" + name + "=2147483648"}), po::error);
            EXPECT_THROW(parse({"--" + name + "=-2147483649"}), po::error);
        };
        (check(option), ...);
    }, AlignOptions::starOptions);
    EXPECT_NO_THROW(AlignParameters(parse({"--aligner=star", "--star_max_multimaps=1",
        "--star_min_junction_overhang=1", "--star_max_segment_gap=0", "--star_min_nonchimeric_score_drop=0",
        "--star_max_chimeric_score_drop=0", "--star_max_intron_length=0", "--min_fragment_length=1"})));
    EXPECT_NO_THROW(AlignParameters(parse({"--aligner=star", "--min_fragment_length=2147483647"})));
    EXPECT_THROW(AlignParameters(parse({"--aligner=star", "--min_fragment_length=0"})), std::invalid_argument);
    EXPECT_THROW(AlignParameters(parse({"--aligner=star", "--min_fragment_length=2147483648"})), std::invalid_argument);
    // Preserve segemehl's existing size_t bounds, even when STAR could not represent them.
    EXPECT_NO_THROW(AlignParameters(parse({"--min_fragment_length=0"})));
    EXPECT_NO_THROW(AlignParameters(parse({"--min_fragment_length=2147483648"})));
    EXPECT_THROW(AlignParameters(parse({"--star_min_junction_overhang=0"})), std::invalid_argument);
    EXPECT_THROW(AlignParameters(parse({"--star_min_junction_overhang=-1"})), std::invalid_argument);
}

TEST_F(ParameterOptionsTests, BackendNamesAreExact) {
    for (const auto& name : {"STAR", "Segemehl", "stars", "", "bwa"}) {
        EXPECT_THROW(AlignParameters(parse({std::string("--aligner=") + name})), std::exception);
    }
}

TEST_F(ParameterOptionsTests, ExplicitStarControlsAreTrackedFromCliAndConfig) {
    std::apply([&](const auto&... option) {
        const auto check = [&](const auto& value) {
            const auto name = value.getLongName();
            SCOPED_TRACE(name);
            AlignParameters cli{parse({"--" + name + "=10"})};
            EXPECT_THROW(cli.validateBackendAvailability(), std::invalid_argument);
            AlignParameters config{parse({}, name + " = 10\n")};
            EXPECT_THROW(config.validateBackendAvailability(), std::invalid_argument);
        };
        (check(option), ...);
    }, AlignOptions::starOptions);
    AlignParameters defaults{parse({}, "aligner = segemehl\n")};
    EXPECT_NO_THROW(defaults.validateBackendAvailability());
}

TEST_F(ParameterOptionsTests, CliOverridesConfigAndConfigOverridesDefaults) {
    AlignParameters p{parse({"--star_max_multimaps=13", "--allow_multimapping=false"},
        "aligner = star\nstar_max_multimaps = 7\nallow_multimapping = true\nstar_max_segment_gap = 6\n")};
    EXPECT_EQ(p.aligner, AlignmentBackend::Star);
    EXPECT_EQ(p.starMaxMultimaps, 13);
    EXPECT_EQ(p.starMaxSegmentGap, 6);
    EXPECT_FALSE(p.multimapAlignments);
    EXPECT_EQ(p.effectiveStarMaxMultimaps(), 1);
}

TEST_F(ParameterOptionsTests, ExplicitBooleansAndInverseFlagsKeepTheirSemantics) {
    const auto inverseOptions = std::make_tuple(
        std::pair{PreprocessOptions::enablePreprocess, "no_preprocess"},
        std::pair{PreprocessOptions::enableDeduplicate, "no_deduplicate"},
        std::pair{PreprocessOptions::trimPolyG, "no_trim_poly_g"},
        std::pair{GeneralOptions::maskMultiCopyGenes, "no_mask_multicopy_genes"},
        std::pair{AlignOptions::allowMultimap, "no_multimapping"},
        std::pair{DetectOptions::removeAltSplicing, "keep_alt_splicing"});
    std::apply([&](const auto&... pair) {
        const auto check = [&](const auto& value) {
            const auto& option = value.first;
            const auto name = option.getLongName();
            const std::string inverse = value.second;
            SCOPED_TRACE(name);
            EXPECT_TRUE(option.extractValue(parse()));
            EXPECT_TRUE(option.extractValue(parse({"--" + name})));
            EXPECT_FALSE(option.extractValue(parse({"--" + name + "=false"})));
            EXPECT_TRUE(option.extractValue(parse({"--" + name + "=true"}, name + " = false\n")));
            EXPECT_FALSE(option.extractValue(parse({}, name + " = false\n")));
            EXPECT_FALSE(option.extractValue(parse({"--" + inverse}, name + " = true\n")));
            EXPECT_FALSE(option.extractValue(parse({}, inverse + " = true\n")));
            // Inverse true wins, including when it came from config (existing behavior).
            EXPECT_FALSE(option.extractValue(parse({"--" + name + "=true"}, inverse + " = true\n")));
            EXPECT_TRUE(option.extractValue(parse({}, inverse + " = false\n")));
        };
        (check(pair), ...);
    }, inverseOptions);
    for (const auto& option : {DetectOptions::excludeSoftClipping, DetectOptions::filterSplicing,
                               DetectOptions::includeWobble}) {
        EXPECT_FALSE(option.extractValue(parse()));
        EXPECT_TRUE(option.extractValue(parse({"--" + option.getLongName() + "=true"})));
        EXPECT_FALSE(option.extractValue(parse({"--" + option.getLongName() + "=false"},
                                              option.getLongName() + " = true\n")));
    }
}

TEST_F(ParameterOptionsTests, ShortcutsExtractTheSameValues) {
    const auto vm = parse({"-T", root.string(), "-C", root.string(), "-t", "7", "-o", "new output",
                           "-r", "new reference", "-f", (root / "features.gff").string(),
                           "-a", "star", "-c", "my config", "-q", "31", "-l", "27", "-h", "-v"});
    const AlignParameters p{vm};
    EXPECT_EQ(p.treatmentsDir, root);
    EXPECT_EQ(p.controlDir, root);
    EXPECT_EQ(p.threadCount, 7);
    EXPECT_EQ(p.outputDir, "new output");
    EXPECT_EQ(p.referenceGenome, "new reference");
    EXPECT_EQ(p.featuresInPath, root / "features.gff");
    EXPECT_EQ(p.aligner, AlignmentBackend::Star);
    EXPECT_EQ(OtherOptions::configFile.extractValue(vm), "my config");
    EXPECT_EQ(PreprocessOptions::minQual.extractValue(vm), 31);
    EXPECT_EQ(PreprocessOptions::minLen.extractValue(vm), 27);
    EXPECT_TRUE(OtherOptions::printHelp.extractValue(vm));
    EXPECT_TRUE(OtherOptions::printVersion.extractValue(vm));
}
TEST_F(ParameterOptionsTests, EveryCanonicalOptionExtractsEquallyFromCliAndConfig) {
    const auto check = [&](const auto& option, const std::string& input, const auto& expected) {
        SCOPED_TRACE(option.getLongName());
        const auto cli = parse({"--" + option.getLongName() + "=" + input});
        const auto config = parse({}, option.getLongName() + " = " + input + "\n");
        EXPECT_EQ(option.extractValue(cli), expected);
        EXPECT_EQ(option.extractValue(config), expected);
        EXPECT_FALSE(cli.at(option.getLongName()).defaulted());
        EXPECT_FALSE(config.at(option.getLongName()).defaulted());
    };
    check(GeneralOptions::trtms, root.string(), root);
    check(GeneralOptions::ctrls, root.string(), root);
    check(GeneralOptions::out, "new output", "new output");
    check(GeneralOptions::logLevel, "warning", LogLevel::WARNING);
    check(GeneralOptions::threads, "7", 7);
    check(GeneralOptions::featuresPath, (root / "features.gff").string(), root / "features.gff");
    check(GeneralOptions::featureTypes, "gene,exon", "gene,exon");
    check(GeneralOptions::featureOrientation, "same", dataTypes::GenomicOrientation::Value::SAME);
    check(GeneralOptions::maskMultiCopyGenes, "false", false);
    check(GeneralOptions::minMultiCopyIdentity, "0.8", 0.8);
    check(GeneralOptions::chunkSize, "13", 13);
    check(PreprocessOptions::enablePreprocess, "false", false);
    check(PreprocessOptions::enableDeduplicate, "false", false);
    check(PreprocessOptions::trimPolyG, "false", false);
    check(PreprocessOptions::minPolyGCount, "6", 6);
    check(PreprocessOptions::adpt5f, "ACGA", "ACGA");
    check(PreprocessOptions::adpt5r, "ACGC", "ACGC");
    check(PreprocessOptions::adpt3f, "ACGG", "ACGG");
    check(PreprocessOptions::adpt3r, "ACGT", "ACGT");
    check(PreprocessOptions::mtrim, "0.2", 0.2);
    check(PreprocessOptions::minOvlTrim, "7", 7);
    check(PreprocessOptions::minQual, "31", 31);
    check(PreprocessOptions::minLen, "27", 27);
    check(PreprocessOptions::wqual, "25", 25);
    check(PreprocessOptions::wtrim, "3", 3);
    check(PreprocessOptions::minOvl, "9", 9);
    check(PreprocessOptions::mmerge, "0.3", 0.3);
    check(AlignOptions::alignmentIndex, "index with spaces",
          std::optional<std::filesystem::path>{"index with spaces"});
    check(AlignOptions::refGenome, "reference with spaces.fa", "reference with spaces.fa");
    check(AlignOptions::aligner, "star", "star");
    check(AlignOptions::allowMultimap, "false", false);
    check(AlignOptions::accuracy, "91", 91);
    check(AlignOptions::minFragmentScore, "19", 19);
    check(AlignOptions::minAlignLength, "21", 21);
    check(AlignOptions::minFragmentLength, "11", 11);
    check(AlignOptions::minSpliceCoverage, "81", 81);
    check(AlignOptions::starMaxMultimaps, "13", 13);
    check(AlignOptions::starMinJunctionOverhang, "14", 14);
    check(AlignOptions::starMaxSegmentGap, "7", 7);
    check(AlignOptions::starMinNonchimericScoreDrop, "11", 11);
    check(AlignOptions::starMaxChimericScoreDrop, "31", 31);
    check(AlignOptions::starMaxIntronLength, "1234", 1234);
    check(DetectOptions::maxPrimaryAlignmentCount, "6", 6);
    check(DetectOptions::minMappingQuality, "9", 9);
    check(DetectOptions::minComplementarity, "0.6", 0.6);
    check(DetectOptions::siteLengthRatio, "0.2", 0.2);
    check(DetectOptions::minDetectLength, "23", 23);
    check(DetectOptions::maxEnergy, "-1.5", -1.5);
    check(DetectOptions::excludeSoftClipping, "true", true);
    check(DetectOptions::filterSplicing, "true", true);
    check(DetectOptions::removeAltSplicing, "false", false);
    check(DetectOptions::splicingTolerance, "6", 6);
    check(DetectOptions::includeWobble, "true", true);
    check(DetectOptions::minHitGroupContribution, "0.2", 0.2);
    check(AnalyzeOptions::maxSelfOverlap, "0.3", 0.3);
    check(AnalyzeOptions::clusteringStrandSpecificity, "specific", dataTypes::GenomicStrandSpecificity::SPECIFIC);
    check(AnalyzeOptions::clusteringDistanceTolerance, "-3", -3);
    check(AnalyzeOptions::clusterFractionOverlap, "0.2", 0.2f);
    check(AnalyzeOptions::maxPadjValue, "0.05", 0.05);
    check(AnalyzeOptions::minimumClusterTranscriptContribution, "2.5", 2.5f);
    check(AnalyzeOptions::minimumSupportPerEffectiveBp, "0.02", 0.02);
    check(AnalyzeOptions::maximumCoverageComponents, "2", 2);
    check(AnalyzeOptions::minimumArmBalance, "0.25", 0.25);
    check(PostprocessOptions::minSegmentFractionOverlap, "0.4", 0.4f);
}
}  // namespace

TEST_F(ParameterOptionsTests, AlignmentIndexExtractionAndPrecedence) {
    EXPECT_FALSE(AlignParameters{parse()}.alignmentIndex);
    for (const auto& flag : {"--alignment_index", "-i"}) {
        const AlignParameters p{parse({flag, "index with spaces", "--no_mask_multicopy_genes"},
                                      "alignment_index = config index\n")};
        EXPECT_EQ(p.alignmentIndex, std::filesystem::path{"index with spaces"});
    }
    const AlignParameters p{parse({}, "alignment_index = config index\nmask_multicopy_genes = false\n")};
    EXPECT_EQ(p.alignmentIndex, std::filesystem::path{"config index"});
    EXPECT_THROW(AlignParameters{parse({"-i", "index"})}, std::invalid_argument);
    EXPECT_THROW(AlignParameters{parse({"-i", "", "--no_mask_multicopy_genes"})}, std::invalid_argument);
}

TEST_F(ParameterOptionsTests, AlignmentIndexValidationOnlyChecksPathAccessibilityAndType) {
    const auto file = root / "arbitrary index name";
    std::ofstream{file} << "Not a parsed index";
    const auto directory = root / "external STAR directory";
    std::filesystem::create_directory(directory);
    for (const auto& backend : {"segemehl", "star"}) {
        const auto valid = std::string{backend} == "star" ? directory : file;
        const auto wrong = std::string{backend} == "star" ? file : directory;
        const auto params = [&](const auto& path) {
            return AlignParameters{parse({"--aligner", backend, "-i", path.string(),
                                           "--no_mask_multicopy_genes"})};
        };
        EXPECT_NO_THROW(params(valid).validateBackendAvailability());
        EXPECT_THROW(params(wrong).validateBackendAvailability(), std::invalid_argument);
        EXPECT_THROW(params(root / "missing").validateBackendAvailability(), std::invalid_argument);
        const auto link = root / (std::string{backend} + " symlink");
        std::filesystem::create_symlink(valid, link);
        EXPECT_NO_THROW(params(link).validateBackendAvailability());
    }
}
