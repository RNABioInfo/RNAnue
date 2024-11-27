#include "ParameterOptions.hpp"

#include <boost/program_options/options_description.hpp>
#include <cstddef>

#include "Constants.hpp"
#include "Orientation.hpp"

namespace pi = constants::pipelines;

auto ParameterOptions::getSubcallOptions() -> po::options_description {
    po::options_description subcall("Subcall");
    subcall.add_options()("subcall", po::value<std::string>(), pi::SUBCALL_DESCRIPTION.c_str());

    return subcall;
}

auto ParameterOptions::getGeneralOptions() -> po::options_description {
    po::options_description general(pi::GENERAL_DESCRIPTION);
    general.add_options()("trtms,t", po::value<std::string>(),
                          "folder containing the raw reads of the treatments including replicates "
                          "located within sub-folders (condition) (required)");
    general.add_options()("ctrls,s", po::value<std::string>(),
                          "folder containing the the raw reads of the controls including "
                          "replicates located within sub-folders (condition)");
    general.add_options()("outdir,o", po::value<std::string>(),
                          "(output) folder in which the results are stored (required)");
    general.add_options()("loglevel", po::value<std::string>()->default_value("info"),
                          "log level [info, warning, error] (default: info)");
    general.add_options()("threads,p", po::value<int>()->default_value(1),
                          "max number of threads to be used (default: 1)");
    general.add_options()("features,f", po::value<std::string>(),
                          "annotation/features in GFF/GTF format (required)");
    general.add_options()(
        "featuretypes", po::value<std::string>()->default_value(std::string{"transcript"}),
        "feature types to be considered for the analysis, can be specified "
        "as --featuretypes 'gene,rRNA' comma seperated values (default: transcript)");
    general.add_options()(
        "orientation",
        po::value<annotation::Orientation>()->default_value(annotation::Orientation::BOTH),
        "orientation of the features to consider in relation to reads [same, "
        "opposite, both] (default: both)");
    general.add_options()("chunksize", po::value<int>()->default_value(pi::defaultChunkSize),
                          "number of reads processed per chunk in parallel (default: 100000)");

    return general;
}

auto ParameterOptions::getPreprocessOptions() -> po::options_description {
    po::options_description preprocess("Preprocess Pipeline");
    preprocess.add_options()(pi::PREPROCESS.c_str(), po::bool_switch()->default_value(true),
                             "whether to include preprocessing of the raw reads in the "
                             "workflow of RNAnue (default: true)");
    preprocess.add_options()("deduplicate", po::bool_switch()->default_value(pi::deduplicate),
                             "remove duplicate reads based on the sequence (default: true)");
    preprocess.add_options()(
        "trimpolyg", po::bool_switch()->default_value(pi::trimpolyG),
        "trim high quality polyG tails from the reads. Applicable for Illumina "
        "NextSeq reads. (default: false)");
    preprocess.add_options()(
        "adpt5f", po::value<std::string>()->default_value(""),
        "single sequence or file [.fasta] of the adapter sequences to be removed "
        "from the 5' end (forward read if PE)");
    preprocess.add_options()(
        "adpt5r", po::value<std::string>()->default_value(""),
        "single sequence or file [.fasta] of the adapter sequences to be removed "
        "from the 5' end of the reverse read (PE only)");
    preprocess.add_options()(
        "adpt3f", po::value<std::string>()->default_value(""),
        "single sequence or file [.fasta] of the adapter sequences to be removed "
        "from the 3' end (forward read if PE)");
    preprocess.add_options()(
        "adpt3r", po::value<std::string>()->default_value(""),
        "single sequence or file [.fasta] of the adapter sequences to be removed "
        "from the 3' end of the reverse read (PE only)");
    preprocess.add_options()(
        "mtrim", po::value<double>()->default_value(pi::defaultAdapterTrimMissmatchRate, "0.05"),
        "rate of mismatches allowed when aligning adapters to sequences (default: 0.05)");
    preprocess.add_options()("minovltrim",
                             po::value<size_t>()->default_value(pi::defaultAdapterTrimMinOverlap),
                             "minimum length of overlap between adapter and read (default: 5)");
    preprocess.add_options()(
        "minqual,q", po::value<size_t>()->default_value(pi::defaultMinMeanPhreadQuality),
        "lower limit for the mean quality (Phred Quality Score) of the reads (default: 20)");
    preprocess.add_options()("minlen,l",
                             po::value<size_t>()->default_value(pi::defaultMinReadLength),
                             "minimum length of the reads (default: 15)");
    preprocess.add_options()(
        "wqual", po::value<size_t>()->default_value(pi::defaultMinWindowPhredQuality),
        "minimum mean quality for each window (Phred Quality Score) (default: 20)");
    preprocess.add_options()("wtrim", po::value<size_t>()->default_value(pi::defaultWindowTrimSize),
                             "window size for quality trimming from 3' end. Selecting '0' will not "
                             "apply quality trimming (default: 0)");
    preprocess.add_options()("minovl",
                             po::value<size_t>()->default_value(pi::defaultMinOverlapMergeSize),
                             "minimal overlap to merge paired-end reads (default: 5)");
    preprocess.add_options()(
        "mmerge",
        po::value<double>()->default_value(pi::defaultMinOverlapMergeMissmatchRate, "0.05"),
        "rate of mismatches allowed when merging paired end reads (default: 0.05)");

    return preprocess;
}

auto ParameterOptions::getAlignOptions() -> po::options_description {
    po::options_description align("Align Pipeline");
    align.add_options()("dbref", po::value<std::string>(), "reference genome (.fasta) (required)");
    align.add_options()("accuracy", po::value<size_t>()->default_value(pi::defaultAlignAccuracy),
                        "minimum percentage of read matches (default: 90, range: 0-100)");
    align.add_options()("minfragsco",
                        po::value<size_t>()->default_value(pi::defaultMinFragmentScore),
                        "minimum score of a spliced fragment (default: 18)");
    align.add_options()("minfraglen",
                        po::value<size_t>()->default_value(pi::defaultMinFragmentLength),
                        "minimum length of a spliced fragment (default: 20)");
    align.add_options()("minsplicecov",
                        po::value<size_t>()->default_value(pi::defaultMinSpliceCoverage),
                        "minimum coverage for spliced transcripts (default: 80, range:0-100)");

    return align;
}

auto ParameterOptions::getDetectOptions() -> po::options_description {
    po::options_description detect("Detect Pipeline");
    detect.add_options()("mapqmin", po::value<size_t>()->default_value(pi::defaultMinMapq),
                         "minimum quality of the alignments (default: 0)");
    detect.add_options()("cmplmin",
                         po::value<double>()->default_value(pi::defaultMinComplementarity),
                         "complementarity cutoff for split reads (default: 0.0; range: 0.0-1.0)");
    detect.add_options()("sitelenratio",
                         po::value<double>()->default_value(pi::defaultMinSitelenRatio, "0.01"),
                         "aligned portion of the read (default: 0.1, range: 0.0-1.0)");
    detect.add_options()(
        "nrgmax", po::value<double>()->default_value(pi::defaultHybridizationEnergyCutoff),
        "hybridization energy cutoff for split reads (default: 0.0, range: >=0.0)");
    detect.add_options()("exclclipping", po::bool_switch()->default_value(false),
                         "exclude soft clipping from the alignments (default: false)");
    detect.add_options()("splicing", po::bool_switch()->default_value(false),
                         "splicing events are removed in the detection of split reads");
    detect.add_options()("splicingtolerance",
                         po::value<int>()->default_value(pi::defaultSplicingTolerance),
                         "tolerance for splicing events (default: 5)");

    return detect;
}

auto ParameterOptions::getAnalyzeOptions() -> po::options_description {
    po::options_description analysis("Analyze Pipeline");
    analysis.add_options()("maxoverlap",
                           po::value<double>()->default_value(pi::defaultMaxOverlap, "0.05"),
                           "maximum fractional overlap between two clusters (default: 0.5)");
    analysis.add_options()("clustdist",
                           po::value<int>()->default_value(pi::defaultClusterTolerance),
                           "threshold distance at which two clusters are merged into a single "
                           "combined cluster, default is to only merge overlapping or blunt "
                           "ended clusters (default: 0)");
    analysis.add_options()(
        "padj", po::value<double>()->default_value(pi::defaultPAdjCutOff),
        "adjusted p-value threshold for outputting an interaction (default: 1.0)");
    analysis.add_options()("mincount",
                           po::value<size_t>()->default_value(pi::defaultMinClusterCount),
                           "minimum number of reads assigned to an interaction (default: 1)");

    return analysis;
}

auto ParameterOptions::getOtherOptions() -> po::options_description {
    po::options_description other("Other");
    other.add_options()("version,v", "display the version number");
    other.add_options()("help,h", "display this help message");
    other.add_options()("config,c", po::value<std::string>(),
                        "configuration file that contains the parameters");

    return other;
}
