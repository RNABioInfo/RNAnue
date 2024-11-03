#include <gtest/gtest.h>

// Standard
#include <filesystem>
#include <vector>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>

// Internal
#include "Deduplicator.hpp"
#include "FastqRecord.hpp"

namespace fs = std::filesystem;
using namespace dataTypes;
using namespace pipelines::preprocess;
using namespace seqan3::literals;

auto testDeduplicatorFastqPath() -> std::string {
    return (std::filesystem::path{__FILE__}.parent_path() / "test_data/duplicated_records.fastq")
        .string();
}

struct DeduplicatorTestParams {
    fs::path recordsPath;
    std::vector<FastqRecord> expectedRecords;
};

class DeduplicatorTests : public ::testing::TestWithParam<DeduplicatorTestParams> {};

TEST_P(DeduplicatorTests, DeduplicateBySequence) {
    const auto& params = GetParam();

    Deduplicator deduplicator(Deduplicator::BySequenceConfig{params.recordsPath});
    const auto results = deduplicator.deduplicate();

    ASSERT_EQ(results.size(), params.expectedRecords.size());

    for (const auto& record : results) {
        std::cout << record.id() << "\n";

        bool found = false;

        for (const auto& expectedRecord : params.expectedRecords) {
            if (record.id() == expectedRecord.id()) {
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

const DeduplicatorTestParams deduplicatorTestParams1{.recordsPath = testDeduplicatorFastqPath(),
                                                     .expectedRecords = expectedRecords1};

INSTANTIATE_TEST_SUITE_P(Default, DeduplicatorTests, testing::Values(deduplicatorTestParams1));
