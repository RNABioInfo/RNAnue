#include <gtest/gtest.h>

#include <atomic>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "pipelines/detect/AsyncSplitReadGroupBuffer.hpp"
#include "pipelines/detect/ComplementarityEvaluationStep.hpp"
#include "pipelines/detect/HybridizationEvaluationStep.hpp"
#include "pipelines/detect/SegemehlReadGroupPreprocessor.hpp"
#include "pipelines/preprocess/CheckedFastqReader.hpp"
#include "pipelines/preprocess/DeduplicatorBySequencePairedEnd.hpp"
#include "pipelines/preprocess/DeduplicatorBySequenceSingleEnd.hpp"
#include "utility/ConcurrentInput.hpp"

using namespace dataTypes;
using namespace pipelines::preprocess;
using namespace pipelines::detect;

namespace {
struct TempFiles {
    std::filesystem::path dir =
        std::filesystem::temp_directory_path() / ("RNAnue-correctness-" + helper::getUUID());
    TempFiles() { std::filesystem::create_directories(dir); }
    ~TempFiles() { std::filesystem::remove_all(dir); }
    auto write(std::string name, std::string data) -> std::filesystem::path {
        auto path = dir / name;
        std::ofstream(path) << data;
        return path;
    }
};
auto fq(std::string id, std::string sequence = "ACGT", char quality = 'I') -> std::string {
    return "@" + id + "\n" + sequence + "\n+\n" + std::string(sequence.size(), quality) + "\n";
}
auto sam(std::string cigar, std::string sequence, int xj = 1, int xh = 1, int flag = 0,
         std::string id = "r", int pos = 1) -> std::string {
    return id + "\t" + std::to_string(flag) + "\tchr1\t" + std::to_string(pos) + "\t60\t" + cigar +
           "\t*\t0\t0\t" + sequence + "\t" + std::string(sequence.size(), 'I') +
           "\tHI:i:1\tNM:i:0\tXJ:i:" + std::to_string(xj) + "\tXH:i:" + std::to_string(xh) +
           "\tXX:i:1\tXY:i:" + std::to_string(sequence.size()) + "\n";
}
auto records(std::string body) -> ReadGroup {
    std::istringstream stream("@SQ\tSN:chr1\tLN:100000\n" + body);
    seqan3::sam_file_input input{stream, seqan3::format_sam{}, SamFieldIDs{}};
    ReadGroup result;
    for (auto& record : input) result.push_back(record);
    return result;
}
}  // namespace

TEST(PairIntegrity, CommonNamesAndOrientation) {
    EXPECT_NO_THROW(validateMateNames("id/1", "id/2"));
    EXPECT_NO_THROW(validateMateNames("id 1:N:0:ACGT", "id 2:N:0:ACGT"));
    EXPECT_NO_THROW(validateMateNames("id comment", "id other comment"));
    EXPECT_THROW(validateMateNames("id/2", "id/1"), std::runtime_error);
    EXPECT_THROW(validateMateNames("id/1 2:N:0:A", "id/2"), std::runtime_error);
    EXPECT_THROW(validateMateNames("run.1", "run.2"), std::runtime_error);
    EXPECT_THROW(validateMateNames("id_1", "id_2"), std::runtime_error);
}

TEST(PairIntegrity, CountsAndIdsCheckedWithOrWithoutDeduplication) {
    TempFiles tmp;
    auto first = tmp.write("R1.fastq", fq("one/1") + fq("two/1"));
    auto second = tmp.write("R2.fastq", fq("one/2"));
    CheckedFastqPairReader reader(first, second);
    ASSERT_TRUE(reader.next());
    try {
        reader.next();
        FAIL();
    } catch (const std::runtime_error& error) {
        std::string message = error.what();
        EXPECT_NE(message.find("pair 2"), std::string::npos);
        EXPECT_NE(message.find(first.string()), std::string::npos);
        EXPECT_NE(message.find(second.string()), std::string::npos);
        EXPECT_NE(message.find("<EOF>"), std::string::npos);
    }
    EXPECT_THROW(DeduplicatorBySequencePairedEnd::deduplicate(first, second), std::runtime_error);
    tmp.write("R2.fastq", fq("wrong/2") + fq("two/2"));
    EXPECT_THROW(DeduplicatorBySequencePairedEnd::deduplicate(first, second), std::runtime_error);
}

TEST(DeduplicationIdentity, SelectsOrdinalsAndPreservesQualityTieRule) {
    TempFiles tmp;
    auto first = tmp.write("R1.fastq", fq("same", "ACGT", '!') + fq("same", "ACGT") +
                                           fq("same", "ACGT") + fq("same", "AAAA"));
    EXPECT_EQ(DeduplicatorBySequenceSingleEnd::deduplicate(first).validRecordOrdinals,
              (std::set<size_t>{1, 3}));
    auto second = tmp.write("R2.fastq", fq("same", "TGCA", '!') + fq("same", "TGCA") +
                                            fq("same", "TGCA") + fq("same", "TTTT"));
    EXPECT_EQ(DeduplicatorBySequencePairedEnd::deduplicate(first, second).validRecordOrdinals,
              (std::set<size_t>{1, 3}));
}

TEST(ConcurrentInput, ProducerAndConsumerFailuresPropagate) {
    struct Count {
        int n{};
        void operator+=(const Count& other) { n += other.n; }
    };
    for (size_t workers : {1, 3}) {
        utility::ConcurrentInput<int> input(1, [n = 0]() mutable -> std::optional<int> {
            if (++n == 4) throw std::runtime_error("producer failure");
            return n;
        });
        EXPECT_THROW(utility::consumeConcurrently(input, workers,
                                                  [&] {
                                                      Count result;
                                                      for (int value : input) result.n += value;
                                                      return result;
                                                  }),
                     std::runtime_error);
        utility::ConcurrentInput<int> endless(1, []() -> std::optional<int> { return 1; });
        EXPECT_THROW(
            utility::consumeConcurrently(endless, workers,
                                         [&]() -> Count {
                                             for (int value : endless) {
                                                 (void)value;
                                                 throw std::runtime_error("consumer failure");
                                             }
                                             return {};
                                         }),
            std::runtime_error);
    }
}

TEST(ConcurrentInput, GroupsSamAndPropagatesMalformedInput) {
    std::istringstream stream("@SQ\tSN:chr1\tLN:100000\n" + sam("4M", "ACGT", 1, 1, 0, "one") +
                              sam("4M", "ACGT", 1, 1, 0, "one") +
                              sam("4M", "ACGT", 1, 1, 0, "two"));
    record_input_t input{stream, seqan3::format_sam{}};
    auto grouped = input | AsyncSplitReadGroupBuffer(1);
    std::vector<size_t> sizes;
    for (auto& group : grouped) sizes.push_back(group.size());
    EXPECT_EQ(sizes, (std::vector<size_t>{2, 1}));
    std::string malformed = sam("4M", "ACGT");
    malformed.replace(malformed.find("\t0\t"), 3, "\tinvalid\t");
    std::istringstream bad("@SQ\tSN:chr1\tLN:100000\n" + malformed);
    record_input_t badInput{bad, seqan3::format_sam{}};
    auto badGroups = badInput | AsyncSplitReadGroupBuffer(1);
    EXPECT_THROW(
        {
            for (auto& group : badGroups) (void)group;
        },
        std::exception);
}

TEST(CigarCorrectness, CountsReferenceAlignedQueryBases) {
    for (const auto& [cigar, length, aligned, clipped] :
         std::vector<std::tuple<std::string, int, int, int>>{
             {"20M5S", 25, 20, 5},
             {"5S20M10S", 35, 20, 15},
             {"5H5S10=2I3X4D90N7M10S5H", 37, 20, 15},
             {"5S", 5, 0, 5}}) {
        auto input = records(sam(cigar, std::string(length, 'A')));
        EXPECT_EQ(alignmentLength(input.front()), aligned) << cigar;
        EXPECT_EQ(softClippedBaseCount(input.front()), clipped) << cigar;
    }
    SamRecord empty;
    EXPECT_EQ(softClippedBaseCount(empty), 0);
    EXPECT_THROW(alignmentLength(empty), std::invalid_argument);
}

TEST(CigarCorrectness, IndependentReferenceCursorForMultipleSplits) {
    for (bool exclude : {false, true}) {
        auto input = records(sam("5H5S10M2I3D90N10M90N10M5S5H", std::string(42, 'A'), 3, 3));
        SegemehlReadGroupPreprocessor pre{{{0, 1, exclude}, 100}};
        auto contexts = pre(std::move(input));
        ASSERT_EQ(contexts.size(), 1);
        std::visit(
            [&](const auto& context) {
                std::vector<int> starts;
                for (const auto& record : context.group->getRecordContainer())
                    starts.push_back(*record.reference_position());
                EXPECT_EQ(starts, (std::vector<int>{0, 103, 203}));
                EXPECT_TRUE(context.isFailed());
            },
            contexts.front());
        EXPECT_EQ(pre.getMetrics().failures.at(HitGroupFailureReason::UNSUPPORTED_MULTISEGMENT), 1);
    }
}

TEST(DetectionSupport, UnmergedPairsRemainExplicitlyUnsupported) {
    auto input = records(sam("20M", std::string(20, 'A'), 1, 1, 65) +
                         sam("20M", std::string(20, 'A'), 1, 1, 129));
    SegemehlReadGroupPreprocessor pre{{{0, 1, false}, 100}};
    pre(std::move(input));
    EXPECT_EQ(pre.getMetrics().failures.at(HitGroupFailureReason::UNSUPPORTED_PAIRED), 1);
}

TEST(FoldingCorrectness, PreservesEnergyZeroAcceptanceAndRepeatedOwnership) {
    HybridizationEvaluationStep step{{0, false}};
    for (int i = 0; i < 30; ++i) {
        auto input = records(sam("12M", "GGGGAAAACCCC", 2, 1) +
                             sam("12M", "AAAAAAAAAAAA", 2, 1, 0, "r", 1001));
        step(EvaluationContext{
            std::make_unique<ChimericHitGroup>(ChimericRecords{input[0], input[1]}, 1)});
    }
    EXPECT_EQ(step.getMetrics().getPassedCount(), 30);
    EXPECT_EQ(step.getMetrics().getMFEValues().front(), 0);
}

TEST(ComplementarityCorrectness, AmbiguousBasesAreNotPositiveEvidence) {
    auto input = records(sam("20M", std::string(20, 'N'), 2, 1) +
                         sam("20M", std::string(20, 'N'), 2, 1, 0, "r", 1001));
    ComplementarityEvaluationStep step{{0.5, 0.1}};
    step(EvaluationContext{
        std::make_unique<ChimericHitGroup>(ChimericRecords{input[0], input[1]}, 1)});
    EXPECT_EQ(step.getMetrics().getPassedCount(), 0);
}

TEST(ComplementarityCorrectness, CountsColumnsAcrossTwoGapRuns) {
    seqan3::nucleotide_scoring_scheme<int8_t> scheme{seqan3::match_score{1},
                                                     seqan3::mismatch_score{-4}};
    auto dna = [](const std::string& s) {
        seqan3::dna5_vector out;
        for (char c : s) {
            seqan3::dna5 base{};
            seqan3::assign_char_to(c, base);
            out.push_back(base);
        }
        return out;
    };
    const auto first = dna(std::string(16, 'A') + std::string(16, 'C') + std::string(16, 'G'));
    const auto second =
        dna(std::string(16, 'A') + "TT" + std::string(16, 'C') + "TT" + std::string(16, 'G'));
    const auto results =
        CoOptimalPairwiseAligner{scheme}.getLocalAlignments(std::tie(first, second));
    ASSERT_FALSE(results.empty());
    EXPECT_NEAR(results.front().complementarity, 48.0 / 52.0, 1e-6);
}

TEST(ComplementarityCorrectness, PreservesWobblePairing) {
    auto input = records(sam("20M", std::string(20, 'G'), 2, 1) +
                         sam("20M", std::string(20, 'T'), 2, 1, 0, "r", 1001));
    ComplementarityEvaluationStep step{{0.5, 0.1}};
    step(EvaluationContext{
        std::make_unique<ChimericHitGroup>(ChimericRecords{input[0], input[1]}, 1)});
    EXPECT_EQ(step.getMetrics().getPassedCount(), 1);
}
