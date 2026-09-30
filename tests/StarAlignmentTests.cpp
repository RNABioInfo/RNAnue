#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include "AlignmentProvenance.hpp"
#include "AlignmentFileInput.hpp"
#include "SamFileUtility.hpp"
#include "ExternalProcess.hpp"
#include "StarAlignmentDecoder.hpp"
#include "StarReadGroupPreprocessor.hpp"
#include "ReadGroupPostScoringStep.hpp"
#include "Utility.hpp"

namespace {
using namespace dataTypes;
using namespace pipelines::align;
using namespace pipelines::detect;
auto sam(const std::string& cigar, const std::string& md, int nm = 0, int flags = 0, int pos = 101,
         std::string ref = "chr1", std::string sequence = std::string(60, 'A'),
         std::string extra = "", int hi = 1, int nh = 1) -> std::string {
    return "read\t" + std::to_string(flags) + "\t" + ref + "\t" + std::to_string(pos) + "\t255\t" +
           cigar + "\t*\t0\t0\t" + sequence + "\t" + std::string(sequence.size(), 'I') +
           "\tHI:i:" + std::to_string(hi) + "\tNH:i:" + std::to_string(nh) +
           "\tAS:i:30\tNM:i:" + std::to_string(nm) + "\tMD:Z:" + md + extra + "\n";
}
auto parse(const std::string& body) -> std::vector<SamRecord> {
    std::istringstream stream{"@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n" + body};
    seqan3::sam_file_input input{stream, seqan3::format_sam{}, SamFieldIDs{}};
    std::vector<SamRecord> result;
    for (auto& record : input) result.push_back(record);
    return result;
}
auto chimera(int secondFlags = 2048) -> std::vector<SamRecord> {
    const auto secondCigar = secondFlags & 16 ? "30M30S" : "30S30M";
    const auto strand = secondFlags & 16 ? "-" : "+";
    return parse(
        sam("30M30S", "10A19", 1, 0, 101, "chr1", std::string(60, 'A'),
            "\tch:A:1\tSA:Z:chr2,201," + std::string(strand) + "," + secondCigar + ",255,0;") +
        sam(secondCigar, "30", 0, secondFlags, 201, "chr2",
            std::string(60, secondFlags & 16 ? 'T' : 'A'),
            "\tch:A:1\tSA:Z:chr1,101,+,30M30S,255,1;"));
}
StarAlignmentDecoder decoder{{"chr1", "chr2"}};
}  // namespace

TEST(StarDecoder, SupplementaryArmsDoNotDuplicateSoftClippedSequence) {
    for (bool exclude : {false, true}) {
        const auto groups = decoder.decode(chimera(), exclude);
        ASSERT_EQ(groups.size(), 1);
        ASSERT_EQ(groups[0].fragments.size(), 2);
        EXPECT_EQ(groups[0].coveredBases[0], 60);
        EXPECT_EQ(groups[0].queryLengths[0], 60);
        EXPECT_EQ(groups[0].fragments[0].sequence().size(), 30);
        EXPECT_EQ(groups[0].fragments[1].sequence().size(), 30);
        EXPECT_EQ(groups[0].fragments[0].tags().get<"NM"_tag>(), 1);
        EXPECT_EQ(groups[0].fragments[1].tags().get<"NM"_tag>(), 0);
        EXPECT_FALSE(groups[0].fragments[0].tags().contains("SA"_tag));
        EXPECT_FALSE(groups[0].fragments[0].tags().contains("MD"_tag));
    }
}

TEST(StarDecoder, ReverseArmUsesOriginalQueryCoordinates) {
    const auto groups = decoder.decode(chimera(2064), false);
    ASSERT_EQ(groups[0].fragments.size(), 2);
    EXPECT_EQ(groups[0].coveredBases[0], 60);
    EXPECT_EQ(groups[0].fragments[1].tags().get<"XX"_tag>(), 30);
    EXPECT_EQ(groups[0].fragments[1].tags().get<"XY"_tag>(), 60);
    EXPECT_EQ(groups[0].fragments[1].sequence().size(), 30);
    EXPECT_TRUE(
        static_cast<bool>(groups[0].fragments[1].flag() & seqan3::sam_flag::on_reverse_strand));
}

TEST(StarDecoder, IntronsAndIndelsHaveSegmentSpecificEditDistances) {
    for (int flag : {0, 16}) {
        const auto raw =
            parse(sam("2S4M2I2M3N3M1D3M2S", "2A6^A3", 4, flag, 101, "chr1", std::string(18, 'A')));
        const auto trimmed = decoder.decode(raw, true);
        ASSERT_EQ(trimmed[0].fragments.size(), 2);
        EXPECT_EQ(trimmed[0].coveredBases[0], 14);
        EXPECT_EQ(trimmed[0].fragments[0].sequence().size(), 8);
        EXPECT_EQ(trimmed[0].fragments[1].sequence().size(), 6);
        EXPECT_EQ(trimmed[0].fragments[0].tags().get<"NM"_tag>(), 3);
        EXPECT_EQ(trimmed[0].fragments[1].tags().get<"NM"_tag>(), 1);
        EXPECT_EQ(trimmed[0].fragments[1].reference_position(), 109);
        const auto extended = decoder.decode(raw, false);
        EXPECT_EQ(extended[0].fragments[0].sequence().size(), 10);
        EXPECT_EQ(extended[0].fragments[1].sequence().size(), 8);
        EXPECT_EQ(alignmentLength(extended[0].fragments[0]), 6);
    }
}

TEST(StarDecoder, AlternativesRemainSeparate) {
    const auto raw = parse(sam("60M", "60", 0, 0, 101, "chr1", std::string(60, 'A'), "", 1, 2) +
                           sam("60M", "60", 0, 256, 201, "chr1", std::string(60, 'A'), "", 2, 2));
    const auto groups = decoder.decode(raw, true);
    ASSERT_EQ(groups.size(), 2);
    EXPECT_EQ(groups[0].rawIndices, std::vector<size_t>{0});
    EXPECT_EQ(groups[1].rawIndices, std::vector<size_t>{1});
}

TEST(StarDecoder, MalformedRelationshipsAndTagsThrowControlledErrors) {
    auto raw = chimera();
    raw[0].tags().erase("HI"_tag);
    EXPECT_THROW(decoder.decode(raw, true), std::runtime_error);
    raw = chimera();
    raw.pop_back();
    EXPECT_THROW(decoder.decode(raw, true), std::runtime_error);
    raw = chimera();
    raw[1].tags()["SA"_tag] = std::string("chr1,999,+,30M30S,255,1;");
    EXPECT_THROW(decoder.decode(raw, true), std::runtime_error);
    raw = chimera();
    raw[0].tags()["NM"_tag] = int32_t{2};
    EXPECT_THROW(decoder.decode(raw, true), std::runtime_error);
    raw = parse(sam("60M", "59"));
    EXPECT_THROW(decoder.decode(raw, true), std::runtime_error);
}

TEST(StarDecoder, OverlappingSegmentsDoNotDoubleCountCoverage) {
    const auto raw = parse(sam("40M20S", "40", 0, 0, 101, "chr1", std::string(60, 'A'),
                               "\tch:A:1\tSA:Z:chr2,201,+,30S30M,255,0;") +
                           sam("30S30M", "30", 0, 2048, 201, "chr2", std::string(60, 'A'),
                               "\tch:A:1\tSA:Z:chr1,101,+,40M20S,255,0;"));
    const auto groups = decoder.decode(raw, false);
    EXPECT_EQ(groups[0].coveredBases[0], 60);
    EXPECT_EQ(groups[0].fragments[0].sequence().size(), 40);
    EXPECT_EQ(groups[0].fragments[1].sequence().size(), 30);
}

TEST(StarDecoder, InternalGapIsNotRetainedAsTerminalClipping) {
    const auto raw = parse(sam("28M32S", "28", 0, 0, 101, "chr1", std::string(60, 'A'),
                               "\tch:A:1\tSA:Z:chr2,201,+,31S29M,255,0;") +
                           sam("31S29M", "29", 0, 2048, 201, "chr2", std::string(60, 'A'),
                               "\tch:A:1\tSA:Z:chr1,101,+,28M32S,255,0;"));
    const auto groups = decoder.decode(raw, false);
    EXPECT_EQ(groups[0].coveredBases[0], 57);
    EXPECT_EQ(groups[0].fragments[0].sequence().size(), 28);
    EXPECT_EQ(groups[0].fragments[1].sequence().size(), 29);
}

TEST(StarPreprocessor, UnsupportedEvidenceAndMapqPolicyArePreserved) {
    SegemehlReadGroupPreprocessorConfig config{{0, 1, true}, 5};
    StarReadGroupPreprocessor preprocessor{config, {"chr1", "chr2"}};
    auto pairs = preprocessor(parse(sam("60M", "60", 0, 65) + sam("60M", "60", 0, 129, 201)));
    ASSERT_EQ(pairs.size(), 1);
    EXPECT_TRUE(
        preprocessor.getMetrics().failures.contains(HitGroupFailureReason::UNSUPPORTED_PAIRED));
    auto multi = preprocessor(parse(sam("20M10N20M10N20M", "60")));
    ASSERT_EQ(multi.size(), 1);
    EXPECT_TRUE(preprocessor.getMetrics().failures.contains(
        HitGroupFailureReason::UNSUPPORTED_MULTISEGMENT));
    config.maxPrimaryAlignmentCount = 1;
    StarReadGroupPreprocessor capped{config, {"chr1", "chr2"}};
    (void)capped(parse(sam("60M", "60", 0, 0, 101, "chr1", std::string(60, 'A'), "", 1, 2) +
                       sam("60M", "60", 0, 256, 201, "chr1", std::string(60, 'A'), "", 2, 2)));
    EXPECT_EQ(capped.getMetrics().failures.at(HitGroupFailureReason::MULTIMAPPING), 2);
}

TEST(StarProvenance, UnknownAndConflictingBackendsAreRejected) {
    seqan3::sam_file_header<> header;
    EXPECT_EQ(utility::alignmentBackend(header), AlignmentBackend::Segemehl);
    utility::stampAlignment(header, AlignmentBackend::Star, "STAR --version");
    EXPECT_EQ(utility::alignmentBackend(header), AlignmentBackend::Star);
    header.comments.push_back("RNAnue alignment backend=segemehl;adapter=1");
    EXPECT_THROW(utility::alignmentBackend(header), std::runtime_error);
    header.comments = {"RNAnue alignment backend=star;adapter=99;version=2.7.11b"};
    EXPECT_THROW(utility::alignmentBackend(header), std::runtime_error);
}

TEST(ExternalProcess, ArgumentsFailuresAndTemporaryCleanup) {
    std::filesystem::path path;
    {
        utility::TemporaryDirectory workspace{std::filesystem::temp_directory_path()};
        path = workspace.path();
        EXPECT_EQ(utility::runExternalProcess("/bin/echo", {"spaces $HOME ; literal"}, path),
                  "spaces $HOME ; literal\n");
        EXPECT_THROW(utility::runExternalProcess("/bin/false", {}, path), std::runtime_error);
        EXPECT_THROW(utility::runExternalProcess("/no/such/executable", {}, path),
                     std::runtime_error);
    }
    EXPECT_FALSE(std::filesystem::exists(path));
}

TEST(StarProvenance, MergeRemapsReferencesAndPreservesMetadataMatesAndQualities) {
    utility::TemporaryDirectory workspace{std::filesystem::temp_directory_path()};
    const auto first = workspace.path() / "first.sam";
    const auto second = workspace.path() / "second.sam";
    const auto output = workspace.path() / "merged.bam";
    const std::string marker = "@CO\tRNAnue alignment backend=star;adapter=1;version=2.7.11b\n";
    const auto write = [&](const auto& path, bool reorder) {
        std::ofstream stream{path};
        stream << (reorder ? "@SQ\tSN:chr2\tLN:1000\n@SQ\tSN:chr1\tLN:1000\n"
                           : "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n")
               << marker << "@PG\tID:aligner\tPN:STAR\tVN:2.7.11b\n";
        stream << "pair\t65\tchr1\t11\t255\t4M\tchr2\t21\t0\tACGT\t~~~~\tHI:i:1\tNH:i:1\tNM:i:"
                  "0\tMD:Z:4\tAS:i:4\tPG:Z:aligner\n";
    };
    write(first, false);
    write(second, true);
    helper::mergeSamFiles({first, second}, output, std::nullopt);
    SamFileUtility::sortByQueryName(output, output, 1);
    utility::AlignmentFileInput input{output};
    EXPECT_EQ(utility::alignmentBackend(input.header()), AlignmentBackend::Star);
    size_t count = 0;
    for (const auto& record : input) {
        ++count;
        EXPECT_EQ(record.tags().get<"PG"_tag>(), count == 1 ? "aligner" : "aligner.1");
        EXPECT_EQ(record.reference_id(), 0);
        EXPECT_EQ(record.mate_reference_id(), 1);
        EXPECT_EQ(record.mate_position(), 20);
        EXPECT_EQ(seqan3::to_phred(record.base_qualities().front()), 93);
    }
    EXPECT_EQ(count, 2);
}

TEST(StarDecoder, ReferenceBoundsAndChimericFlagsAreValidated) {
    StarAlignmentDecoder bounded{{"chr1", "chr2"}, {150, 1000}};
    EXPECT_THROW(bounded.decode(parse(sam("60M", "60")), true), std::runtime_error);
    auto raw = chimera();
    raw[1].flag() = seqan3::sam_flag::none;
    EXPECT_THROW(decoder.decode(raw, true), std::runtime_error);
    raw = chimera();
    raw[0].tags()["ch"_tag] = '0';
    EXPECT_THROW(decoder.decode(raw, true), std::runtime_error);
}
