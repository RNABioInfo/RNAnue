#include <gtest/gtest.h>

// Standard
#include <filesystem>
#include <string>
#include <vector>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>

// Internal
#include "DeduplicationConfig.hpp"
#include "DeduplicationOutput.hpp"
#include "Deduplicator.hpp"
#include "FastqRecord.hpp"
#include "TestFilePath.hpp"

namespace fs = std::filesystem;
using namespace dataTypes;
using namespace pipelines::preprocess;
using namespace seqan3::literals;

struct DeduplicatorSingleTestParams {
    fs::path recordsPath;
    std::vector<FastqRecord> expectedRecords;
};

class DeduplicatorSingleTests : public ::testing::TestWithParam<DeduplicatorSingleTestParams> {};

TEST_P(DeduplicatorSingleTests, DeduplicateBySequenceSingle) {
    const auto& params = GetParam();

    const auto results =
        Deduplicator::deduplicate(DeduplicationBySequenceSingleConfig{params.recordsPath});

    ASSERT_EQ(results.validRecordIDs.size(), params.expectedRecords.size());

    for (const auto& recordID : results.validRecordIDs) {
        std::cout << recordID << "\n";

        bool found = false;

        for (const auto& expectedRecord : params.expectedRecords) {
            if (recordID == expectedRecord.id()) {
                found = true;
                break;
            }
        }

        EXPECT_TRUE(found);
    }
}

const std::vector<FastqRecord> expectedRecords1{
    {"AGGTGGAGTCGACGTATAAGCCGGGTTCTGTTCCGC"_dna5, "SRR18331301.1 1 length=36",
     "GGIIGIIIGGIGIIIIIGIAGGGGGGGIIGIGIIGI"_phred42},
    {"CAGCCATCTTCAAAACGATGCTGACCTCGTGACCAA"_dna5, "SRR18331301.4 4 length=36",
     "GIIIIGIIIIIIGIIIIIIIIIIIGGIIIIIGGIIJ"_phred42}};

const DeduplicatorSingleTestParams deduplicatorTestParams1{
    .recordsPath = getTestFilePath("duplicatedRecords.fastq"), .expectedRecords = expectedRecords1};

INSTANTIATE_TEST_SUITE_P(Default, DeduplicatorSingleTests,
                         testing::Values(deduplicatorTestParams1));

struct DeduplicatorPairedTestParams {
    fs::path recordsPathFwd;
    fs::path recordsPathRev;
    std::vector<std::pair<FastqRecord, FastqRecord>> expectedRecordPairs;
};

class DeduplicatorPairedTests : public ::testing::TestWithParam<DeduplicatorPairedTestParams> {};

TEST_P(DeduplicatorPairedTests, DeduplicateBySequencePaired) {
    const auto& params = GetParam();

    const auto results = Deduplicator::deduplicate(DeduplicationBySequencePairedConfig{
        .recordsPathFwd = params.recordsPathFwd, .recordsPathRev = params.recordsPathRev});

    ASSERT_EQ(results.validRecordIDs.size(), params.expectedRecordPairs.size());

    for (const std::string& recordID : results.validRecordIDs) {
        bool found = false;

        for (const auto& expectedRecord : params.expectedRecordPairs) {
            if (recordID == expectedRecord.first.id() && recordID == expectedRecord.second.id()) {
                found = true;
                break;
            }
        }

        EXPECT_TRUE(found);
    }
}

const std::vector<std::pair<FastqRecord, FastqRecord>> expectedRecords2{
    {{"AGGTGGAGTCGACGTATAAGCCGGGTTCTGTTCCGC"_dna5, "SRR18331301.1 1 length=36",
      "GGIIGIIIGGIGIIIIIGIAGGGGGGGIIGIGIIGI"_phred42},
     {"AGGTGGAGTCGACGTATAAGCCGGGTTCTGTTCCGC"_dna5, "SRR18331301.1 1 length=36",
      "GGIIGIIIGGIGIIIIIGIAGGGGGGGIIGIGIIGI"_phred42}},
    {{"AGGTGGAGTCGACGTATAAGCCGGGTTCTGTTCCGC"_dna5, "SRR18331301.2 2 length=36",
      "GGIIGIIIGGIGIIIIIGIAGGGGGGGIIGIGIIGI"_phred42},
     {"AGGTGGAGTCGACGTATAAGCCGGGTTCTGTTCCGA"_dna5, "SRR18331301.2 2 length=36",
      "GGIIGIIIGGIGIIIIIGIAGGGGGGGIIGIGIIGI"_phred42}},
    {{"CAGCCATCTTCAAAACGATGCTGACCTCGTGACCAA"_dna5, "SRR18331301.4 4 length=36",
      "GIIIIGIIIIIIGIIIIIIIIIIIGGIIIIIGGIIJ"_phred42},
     {"CAGCCATCTTCAAAACGATGCTGACCTCGTGACCAA"_dna5, "SRR18331301.4 4 length=36",
      "GIIIIGIIIIIIGIIIIIIIIIIIGGIIIIIGGIIJ"_phred42}}};

const DeduplicatorPairedTestParams deduplicatorTestParams2{
    .recordsPathFwd = getTestFilePath("duplicatedRecords.fastq"),
    .recordsPathRev = getTestFilePath("duplicatedRecords2.fastq"),
    .expectedRecordPairs = expectedRecords2};

INSTANTIATE_TEST_SUITE_P(Default, DeduplicatorPairedTests,
                         testing::Values(deduplicatorTestParams2));
