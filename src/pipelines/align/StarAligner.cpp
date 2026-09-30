#include "StarAligner.hpp"

#include <htslib/hts.h>
#include <unistd.h>
#include <zlib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <set>
#include <sstream>
#include <stdexcept>

#include "AlignmentFileInput.hpp"
#include "AlignmentProvenance.hpp"
#include "CheckedFastqReader.hpp"
#include "ExternalProcess.hpp"
#include "Logger.hpp"
#include "SamFileUtility.hpp"
#include "StarAlignmentDecoder.hpp"
#include "StarExecutable.hpp"

namespace pipelines::align {
namespace {
namespace fs = std::filesystem;

auto readText(const fs::path& path) -> std::string {
    std::ifstream input{path, std::ios::binary};
    return {std::istreambuf_iterator<char>{input}, {}};
}

auto fingerprint(const fs::path& path) -> std::string {
    std::ifstream input{path, std::ios::binary};
    if (!input) {
        throw std::runtime_error("Cannot open reference genome: " + path.string());
    }
    std::unique_ptr<hts_md5_context, decltype(&hts_md5_destroy)> md5{hts_md5_init(),
                                                                     hts_md5_destroy};
    if (!md5) {
        throw std::bad_alloc{};
    }

    std::array<char, 65536> buffer;
    while (input.read(buffer.data(), buffer.size()) || (input.gcount() != 0)) {
        hts_md5_update(md5.get(), buffer.data(), static_cast<unsigned long>(input.gcount()));
    }

    if (!input.eof()) {
        throw std::runtime_error("Cannot fingerprint reference: " + path.string());
    }

    unsigned char digest[16];
    char hex[33];

    hts_md5_final(digest, md5.get());
    hts_md5_hex(hex, digest);
    return hex;
}
void decompress(const fs::path& input, const fs::path& output) {
    auto* const handle = gzopen(input.c_str(), "rb");

    if (handle == nullptr) {
        throw std::runtime_error("Cannot open sequence input: " + input.string());
    }

    struct Gzip {
        gzFile handle;
        ~Gzip() { gzclose(handle); }
    }

    owner{handle};

    std::ofstream stream{output, std::ios::binary};
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    std::array<char, 65536> buffer;
    int n;
    while ((n = gzread(handle, buffer.data(), buffer.size())) > 0) {
        stream.write(buffer.data(), n);
    }

    int code{};
    const char* reason = gzerror(handle, &code);

    if (n < 0 || (code != Z_OK && code != Z_STREAM_END)) {
        throw std::runtime_error("Cannot decompress " + input.string() + ": " + reason);
    }
}
// A text reader avoids re-encoding qualities or identifiers while validating STAR input.
struct FastqReader {
    std::ifstream stream;
    fs::path path;
    size_t ordinal{};
    explicit FastqReader(const fs::path& value) : stream(value), path(value) {
        if (!stream) {
            throw std::runtime_error("Cannot open FASTQ: " + value.string());
        }
    }
    auto next() -> std::optional<std::array<std::string, 4>> {
        std::array<std::string, 4> record;
        if (!std::getline(stream, record[0])) {
            if (!stream.eof()) throw std::runtime_error("Cannot read FASTQ: " + path.string());
            return std::nullopt;
        }

        ++ordinal;
        for (size_t i = 1; i < 4; ++i) {
            if (!std::getline(stream, record[i])) {
                throw std::runtime_error("Truncated FASTQ: " + path.string());
            }
        }

        for (auto& line : record) {
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
        }

        if (!record[0].starts_with('@') || record[0].size() < 2 || record[0][1] == ' ' ||
            record[0][1] == '\t' || !record[2].starts_with('+') || record[1].empty() ||
            record[1].size() != record[3].size() ||
            record[1].find_first_not_of("ACGTUNRYKMSWBDHVacgtunrykmswbdhv") != std::string::npos ||
            std::ranges::any_of(record[3], [](unsigned char c) { return c < 33 || c > 126; })) {
            throw std::runtime_error("Malformed FASTQ record " + std::to_string(ordinal) + ": " +
                                     path.string());
        }

        return record;
    }
};
void writeFastq(std::ostream& output, const std::array<std::string, 4>& record) {
    for (const auto& line : record) {
        output << line << '\n';
    }
}
auto commandString(const fs::path& executable, const std::vector<std::string>& args)
    -> std::string {
    std::ostringstream result;
    result << std::quoted(executable.string());
    for (const auto& arg : args) {
        result << ' ' << std::quoted(arg);
    }
    auto text = result.str();
    std::ranges::replace(text, '\t', ' ');
    std::ranges::replace(text, '\n', ' ');
    std::ranges::replace(text, '\r', ' ');
    return text;
}
struct IndexLock {
    fs::path path;
    explicit IndexLock(fs::path value) : path(std::move(value)) {
        if (!fs::create_directory(path)) {
            throw std::runtime_error(
                "STAR index is locked: " + path.string() +
                "; wait for its builder, or remove the lock after verifying no builder is running");
        }
    }
    ~IndexLock() {
        std::error_code e;
        fs::remove(path, e);
    }
};
}  // namespace

void StarAligner::buildIndex() {
    const auto genomePath = fs::absolute(parameters.referenceGenome);
    const auto hash = fingerprint(genomePath);

    reference = dataTypes::SamReference{{}, {}};
    size_t genomeSize = 0;
    std::set<std::string> uniqueIDs;
    seqan3::sequence_file_input genome{genomePath};

    for (const auto& record : genome) {
        const auto id = record.id().substr(0, record.id().find_first_of(" \t"));

        if (id.empty() || !uniqueIDs.insert(id).second || record.sequence().empty()) {
            throw std::runtime_error(
                "STAR reference requires unique identifiers and nonempty sequences");
        }

        reference.referenceIDs.push_back(id);
        reference.referenceLengths.push_back(record.sequence().size());
        genomeSize += record.sequence().size();
    }
    if (genomeSize == 0U) {
        throw std::runtime_error("STAR reference genome is empty");
    }

    const int bases = std::clamp(
        static_cast<int>(std::floor((std::log2(static_cast<double>(genomeSize)) / 2.0) - 1)), 1,
        14);

    const std::string metadata =
        "RNAnue STAR index 1\nversion=" + std::string(utility::star::version()) +
        "\nreference_md5=" + hash + "\ngenomeSAindexNbases=" + std::to_string(bases) +
        "\ngenomeChrBinNbits=18\ngenomeSAsparseD=1\ncomplete=true\n";

    const auto name = genomePath.filename().string() + ".star_index";

    const std::array<fs::path, 2> candidates{genomePath.parent_path() / name,
                                             fs::absolute(parameters.outputDir) / name};

    const auto valid = [&](const fs::path& candidate) {
        if (readText(candidate / "rnanue-index.meta") != metadata) return false;
        for (const auto* const file :
             {"Genome", "SA", "SAindex", "genomeParameters.txt", "chrName.txt", "chrLength.txt",
              "chrStart.txt", "chrNameLength.txt"}) {
            std::error_code error;

            if (!fs::is_regular_file(candidate / file, error) ||
                fs::file_size(candidate / file, error) == 0 || error) {
                return false;
            }
        }
        return true;
    };
    for (const auto& candidate : candidates) {
        if (valid(candidate)) {
            indexPath = candidate;
            Logger::log("Reusing STAR index: ", indexPath);
            return;
        }
    }
    indexPath = access(genomePath.parent_path().c_str(), W_OK) == 0 ? candidates[0] : candidates[1];

    fs::create_directories(indexPath.parent_path());
    IndexLock lock{indexPath.string() + ".lock"};

    if (valid(indexPath)) {
        return;
    }

    utility::TemporaryDirectory workspace{indexPath.parent_path()};
    const auto staged = workspace.path() / "index";

    fs::create_directory(staged);
    auto fasta = genomePath;

    if (genomePath.extension() == ".gz") {
        fasta = workspace.path() / "genome.fa";
        decompress(genomePath, fasta);
    }

    const std::vector<std::string> args{"--runMode",
                                        "genomeGenerate",
                                        "--genomeDir",
                                        staged.string(),
                                        "--genomeFastaFiles",
                                        fasta.string(),
                                        "--genomeSAindexNbases",
                                        std::to_string(bases),
                                        "--runThreadN",
                                        std::to_string(parameters.threadCount),
                                        "--outFileNamePrefix",
                                        (workspace.path() / "index.").string()};

    Logger::log("Building STAR index: ", indexPath);
    (void)utility::runExternalProcess(utility::star::executablePath(), args, workspace.path());
    // STAR also copies its generation log into genomeDir, regardless of the prefix.
    fs::remove(staged / "Log.out");
    if (fingerprint(genomePath) != hash) {
        throw std::runtime_error("Reference changed while STAR was indexing it");
    }

    {
        std::ofstream meta{staged / "rnanue-index.meta"};
        meta.exceptions(std::ios::badbit | std::ios::failbit);
        meta << metadata;
    }

    if (!valid(staged)) {
        throw std::runtime_error("STAR generated an incomplete index");
    }

    const auto previous = workspace.path() / "previous";
    const bool existed = fs::exists(indexPath);
    if (existed) {
        fs::rename(indexPath, previous);
    }

    try {
        fs::rename(staged, indexPath);
    } catch (...) {
        if (existed) {
            fs::rename(previous, indexPath);
        }
        throw;
    }
}

auto StarAligner::alignmentArguments() const -> std::vector<std::string> {
    return {"--genomeDir",
            indexPath.string(),
            "--runThreadN",
            std::to_string(parameters.threadCount),
            "--genomeLoad",
            "NoSharedMemory",
            "--twopassMode",
            "None",
            "--outSAMtype",
            "BAM",
            "Unsorted",
            "--chimOutType",
            "WithinBAM",
            "SoftClip",
            "--outSAMattributes",
            "NH",
            "HI",
            "AS",
            "NM",
            "MD",
            "ch",
            "--readNameSeparator",
            "space",
            "--outSJtype",
            "None",
            "--quantMode",
            "-",
            "--outSAMunmapped",
            "None",
            "--outFilterScoreMinOverLread",
            "0",
            "--outFilterMatchNminOverLread",
            "0",
            "--chimSegmentMin",
            std::to_string(parameters.minimumFragmentLength),
            "--outFilterMultimapNmax",
            std::to_string(parameters.effectiveStarMaxMultimaps()),
            "--chimMultimapNmax",
            std::to_string(parameters.effectiveStarMaxMultimaps()),
            "--chimJunctionOverhangMin",
            std::to_string(parameters.effectiveStarMinJunctionOverhang()),
            "--chimSegmentReadGapMax",
            std::to_string(parameters.starMaxSegmentGap),
            "--chimNonchimScoreDropMin",
            std::to_string(parameters.starMinNonchimericScoreDrop),
            "--chimScoreDropMax",
            std::to_string(parameters.starMaxChimericScoreDrop),
            "--alignIntronMax",
            std::to_string(parameters.starMaxIntronLength)};
}

void StarAligner::alignSingleReads(const fs::path& input, const fs::path& output) const {
    align(input, std::nullopt, output);
}
void StarAligner::alignPairedReads(const fs::path& forward, const fs::path& reverse,
                                   const fs::path& output) const {
    align(forward, reverse, output);
}

void StarAligner::align(const fs::path& forward, const std::optional<fs::path>& reverse,
                        const fs::path& output) const {
    utility::TemporaryDirectory workspace{output.parent_path()};
    const auto rawForward = workspace.path() / "forward.raw.fastq";
    const auto rawReverse = workspace.path() / "reverse.raw.fastq";
    decompress(forward, rawForward);

    if (reverse) {
        decompress(*reverse, rawReverse);
    }

    FastqReader first{rawForward};
    std::optional<FastqReader> second;
    if (reverse) {
        second.emplace(rawReverse);
    }

    // Single reads can be filtered before mapping; paired records are retained together and
    // filtered per mate after mapping, preserving their paired status for detection.
    const auto firstPath = workspace.path() / "forward.fastq";
    const auto secondPath = workspace.path() / "reverse.fastq";
    size_t entries = 0;
    {
        std::ofstream f{firstPath};
        std::ofstream r{secondPath};

        f.exceptions(std::ios::badbit | std::ios::failbit);
        r.exceptions(std::ios::badbit | std::ios::failbit);

        while (true) {
            auto a = first.next();
            auto b = second ? second->next() : std::nullopt;

            if (second && static_cast<bool>(a) != static_cast<bool>(b)) {
                throw std::runtime_error("STAR paired inputs contain unequal record counts");
            }

            if (!a) {
                break;
            }

            if (b) {
                preprocess::validateMateNames((*a)[0].substr(1), (*b)[0].substr(1));
            }

            if ((!b && (*a)[1].size() < parameters.minLengthThreshold) ||
                (b && (*a)[1].size() < parameters.minLengthThreshold &&
                 (*b)[1].size() < parameters.minLengthThreshold)) {
                continue;
            }

            writeFastq(f, *a);

            if (b) {
                writeFastq(r, *b);
            }
            ++entries;
        }
    }
    auto args = alignmentArguments();
    args.insert(args.end(), {"--readFilesIn", firstPath.string()});

    if (reverse) {
        args.push_back(secondPath.string());
    }

    args.insert(args.end(), {"--outFileNamePrefix", (workspace.path() / "mapping.").string()});
    const auto executable = utility::star::executablePath();
    const auto command = commandString(executable, args);
    const auto finalPath = workspace.path() / "final.bam";

    if (entries == 0U) {
        SamFileUtility::writeHeaderOnlyFile(finalPath, reference);
        utility::stampAlignmentFile(finalPath, AlignmentBackend::Star, command);
    } else {
        (void)utility::runExternalProcess(executable, args, workspace.path());
        const auto raw = workspace.path() / "mapping.Aligned.out.bam";
        const auto rawInspection = SamFileUtility::inspect(raw);

        if (!rawInspection.isReadable()) {
            throw std::runtime_error("STAR did not produce a readable BAM");
        }

        SamFileUtility::sortByQueryName(raw, raw, parameters.threadCount);

        {
            utility::AlignmentFileInput input{raw};
            dataTypes::SamReference dictionary{input.header()};
            seqan3::sam_file_output out{finalPath, dictionary.referenceIDs,
                                        dictionary.referenceLengths};
            out.header().program_infos = input.header().program_infos;
            utility::stampAlignment(out.header(), AlignmentBackend::Star, command);
            out.header().sorting = "queryname";
            StarAlignmentDecoder decoder{dictionary.referenceIDs, dictionary.referenceLengths};

            std::vector<dataTypes::SamRecord> group;
            std::vector<typename decltype(input)::record_type> originals;

            const auto flush = [&] {
                for (const auto& hypothesis : decoder.decode(group, true)) {
                    std::array<bool, 2> keep{};

                    for (size_t m = 0; m < 2; ++m) {
                        const auto length = hypothesis.queryLengths[m];
                        // Avoid overflow while preserving segemehl's integer percent rounding.
                        const auto required =
                            ((length / 100) * parameters.minimumSpliceCoverage) +
                            ((length % 100) * parameters.minimumSpliceCoverage / 100);
                        keep[m] = length >= parameters.minLengthThreshold &&
                                  hypothesis.coveredBases[m] >= required;
                    }

                    if (hypothesis.fragments.size() > 1 && !hypothesis.paired) {
                        for (const auto& fragment : hypothesis.fragments) {
                            if (fragment.sequence().size() < parameters.minimumFragmentLength) {
                                keep[0] = false;
                            }
                        }
                    }

                    for (const auto i : hypothesis.rawIndices) {
                        const auto m = (static_cast<unsigned>(group[i].flag()) & 0x80) ? 1 : 0;

                        if (keep[m]) {
                            out.push_back(originals[i]);
                        }
                    }
                }
                group.clear();
                originals.clear();
            };

            if (rawInspection.hasRecords()) {
                for (auto& record : input) {
                    if (!group.empty() && group.front().id() != record.id()) {
                        flush();
                    }

                    std::vector<seqan3::phred42> qualities;
                    for (const auto quality : record.base_qualities()) {
                        qualities.push_back(
                            seqan3::phred42{}.assign_phred(seqan3::to_phred(quality)));
                    }

                    group.emplace_back(record.id(), record.flag(), record.reference_id(),
                                       record.reference_position(), record.mapping_quality(),
                                       record.cigar_sequence(), record.sequence(),
                                       std::move(qualities), record.tags());
                    originals.push_back(std::move(record));
                }
            }

            if (!group.empty()) {
                flush();
            }
        }
    }

    const auto inspection = SamFileUtility::inspect(finalPath);

    if (!inspection.isReadable() || inspection.hasMissingEof()) {
        throw std::runtime_error("Invalid filtered STAR BAM");
    }

    fs::rename(finalPath, output);
}
}  // namespace pipelines::align
