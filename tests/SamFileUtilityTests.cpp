// NOLINTBEGIN(readability-magic-numbers)

#include <gtest/gtest.h>

// Standard
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <thread>
#include <vector>

// htslib
#include <htslib/hts.h>
#include <htslib/sam.h>

// Internal
#include "SamFileUtility.hpp"
#include "Utility.hpp"

namespace {

namespace fs = std::filesystem;

auto testReference() -> dataTypes::SamReference {
    return dataTypes::SamReference{std::deque<std::string>{"ref"}, std::vector<size_t>{100}};
}

auto uniqueTestDir() -> fs::path {
    static std::atomic_uint64_t counter{0};
    const auto timestamp = std::chrono::steady_clock::now().time_since_epoch().count();
    return fs::temp_directory_path() /
           ("rnanue_sam_file_tests_" + std::to_string(timestamp) + "_" +
            std::to_string(counter++));
}

class SamFileUtilityTest : public testing::Test {
   protected:
    SamFileUtilityTest() : tmpDir(uniqueTestDir()) { fs::create_directories(tmpDir); }

    ~SamFileUtilityTest() override {
        std::error_code ignoredError;
        fs::remove_all(tmpDir, ignoredError);
    }

    [[nodiscard]] auto path(const std::string& filename) const -> fs::path {
        return tmpDir / filename;
    }

    fs::path tmpDir;
};

void writeHeaderOnlyBam(const fs::path& path) {
    SamFileUtility::writeHeaderOnlyFile(path, testReference());
}

void writeRecordBam(const fs::path& path, const std::string& qname = "read1") {
    using SamFilePtr = std::unique_ptr<samFile, decltype(&hts_close)>;
    using HeaderPtr = std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)>;
    using RecordPtr = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>;

    SamFilePtr out{sam_open(path.string().c_str(), "wb"), hts_close};
    ASSERT_NE(out, nullptr);

    const std::string headerText = "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref\tLN:100\n";
    HeaderPtr header{sam_hdr_parse(headerText.size(), headerText.c_str()), sam_hdr_destroy};
    ASSERT_NE(header, nullptr);
    ASSERT_EQ(sam_hdr_write(out.get(), header.get()), 0);

    RecordPtr record{bam_init1(), bam_destroy1};
    ASSERT_NE(record, nullptr);

    const uint32_t cigar = bam_cigar_gen(4, BAM_CMATCH);
    ASSERT_GE(bam_set1(record.get(), qname.size() + 1, qname.c_str(), 0, 0, 0, 60, 1,
                       &cigar, -1, -1, 0, 4, "ACGT", "IIII", 0),
              0);
    ASSERT_GE(sam_write1(out.get(), header.get(), record.get()), 0);

    auto* rawOut = out.release();
    ASSERT_EQ(hts_close(rawOut), 0);
}

void removeBgzfEofMarker(const fs::path& path) {
    const auto size = fs::file_size(path);
    ASSERT_GT(size, 28U);
    fs::resize_file(path, size - 28);
}

void truncateMidBlock(const fs::path& path) {
    const auto size = fs::file_size(path);
    ASSERT_GT(size, 48U);
    fs::resize_file(path, size / 2);
}

}  // namespace

TEST_F(SamFileUtilityTest, InspectMissingPath) {
    const auto inspection = SamFileUtility::inspect(path("missing.bam"));

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::Missing);
    EXPECT_FALSE(inspection.isReadable());
}

TEST_F(SamFileUtilityTest, InspectZeroByteFile) {
    std::ofstream{path("zero.bam")};

    const auto inspection = SamFileUtility::inspect(path("zero.bam"));

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::ZeroByte);
    EXPECT_FALSE(inspection.isReadable());
}

TEST_F(SamFileUtilityTest, InspectHeaderOnlyBam) {
    writeHeaderOnlyBam(path("header_only.bam"));

    const auto inspection = SamFileUtility::inspect(path("header_only.bam"));

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::HeaderOnly);
    EXPECT_TRUE(inspection.isReadable());
    EXPECT_FALSE(inspection.hasRecords());
    EXPECT_EQ(inspection.recordCount, 0U);
    EXPECT_EQ(inspection.referenceIDs, std::deque<std::string>{"ref"});
    EXPECT_EQ(inspection.referenceLengths, std::vector<size_t>{100});
}

TEST_F(SamFileUtilityTest, InspectRecordContainingBam) {
    writeRecordBam(path("records.bam"));

    const auto inspection = SamFileUtility::inspect(path("records.bam"));

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::HasRecords);
    EXPECT_TRUE(inspection.isReadable());
    EXPECT_TRUE(inspection.hasRecords());
    EXPECT_EQ(inspection.recordCount, 1U);
}

TEST_F(SamFileUtilityTest, InspectReadableBamMissingBgzfEofMarker) {
    writeRecordBam(path("missing_eof.bam"));
    removeBgzfEofMarker(path("missing_eof.bam"));

    const auto inspection = SamFileUtility::inspect(path("missing_eof.bam"));

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::HasRecordsMissingEof);
    EXPECT_TRUE(inspection.isReadable());
    EXPECT_TRUE(inspection.hasRecords());
    EXPECT_TRUE(inspection.hasMissingEof());
    EXPECT_EQ(inspection.recordCount, 1U);
}

TEST_F(SamFileUtilityTest, InspectTruncatedMidBlockBam) {
    writeRecordBam(path("truncated.bam"));
    truncateMidBlock(path("truncated.bam"));

    const auto inspection = SamFileUtility::inspect(path("truncated.bam"));

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::UnreadableOrTruncated);
    EXPECT_FALSE(inspection.isReadable());
}

TEST_F(SamFileUtilityTest, InspectWithRetriesAcceptsFileThatAppearsLate) {
    const auto delayedPath = path("delayed.bam");
    const auto tempPath = path("delayed.tmp.bam");
    std::thread writer{[delayedPath, tempPath] {
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        SamFileUtility::writeHeaderOnlyFile(tempPath, testReference());
        fs::rename(tempPath, delayedPath);
    }};

    const auto inspection = SamFileUtility::inspectWithRetries(delayedPath, 8, 5);
    writer.join();

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::HeaderOnly);
    EXPECT_TRUE(inspection.isReadable());
}

TEST_F(SamFileUtilityTest, MergeEmptyInputListWithReferenceCreatesHeaderOnlyBam) {
    helper::mergeSamFiles({}, path("merged.bam"), testReference());

    const auto inspection = SamFileUtility::inspect(path("merged.bam"));

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::HeaderOnly);
    EXPECT_EQ(inspection.recordCount, 0U);
    EXPECT_EQ(inspection.referenceIDs, std::deque<std::string>{"ref"});
}

TEST_F(SamFileUtilityTest, MergeHeaderOnlyInputsCreatesHeaderOnlyBam) {
    writeHeaderOnlyBam(path("input_a.bam"));
    writeHeaderOnlyBam(path("input_b.bam"));

    helper::mergeSamFiles({path("input_a.bam"), path("input_b.bam")}, path("merged.bam"),
                          std::nullopt);

    const auto inspection = SamFileUtility::inspect(path("merged.bam"));

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::HeaderOnly);
    EXPECT_EQ(inspection.recordCount, 0U);
    EXPECT_EQ(inspection.referenceIDs, std::deque<std::string>{"ref"});
}

TEST_F(SamFileUtilityTest, MergeMixedHeaderOnlyAndRecordsPreservesRecords) {
    writeHeaderOnlyBam(path("header_only.bam"));
    writeRecordBam(path("records.bam"));

    helper::mergeSamFiles({path("header_only.bam"), path("records.bam")}, path("merged.bam"),
                          std::nullopt);

    const auto inspection = SamFileUtility::inspect(path("merged.bam"));

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::HasRecords);
    EXPECT_EQ(inspection.recordCount, 1U);
    EXPECT_EQ(inspection.referenceIDs, std::deque<std::string>{"ref"});
}

TEST_F(SamFileUtilityTest, MergeReadableMissingEofInputRewritesCleanOutput) {
    writeRecordBam(path("missing_eof.bam"));
    removeBgzfEofMarker(path("missing_eof.bam"));

    helper::mergeSamFiles({path("missing_eof.bam")}, path("merged.bam"), std::nullopt);

    const auto inspection = SamFileUtility::inspect(path("merged.bam"));

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::HasRecords);
    EXPECT_FALSE(inspection.hasMissingEof());
    EXPECT_EQ(inspection.recordCount, 1U);
}

TEST_F(SamFileUtilityTest, SortSameInputOutputKeepsValidBam) {
    writeRecordBam(path("sort_in_place.bam"));

    SamFileUtility::sortByQueryName(path("sort_in_place.bam"), path("sort_in_place.bam"), 1);

    const auto inspection = SamFileUtility::inspect(path("sort_in_place.bam"));

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::HasRecords);
    EXPECT_EQ(inspection.recordCount, 1U);
}

TEST_F(SamFileUtilityTest, SortReadableMissingEofInputRepairsOutput) {
    writeRecordBam(path("sort_missing_eof.bam"));
    removeBgzfEofMarker(path("sort_missing_eof.bam"));

    SamFileUtility::sortByQueryName(path("sort_missing_eof.bam"), path("sort_missing_eof.bam"), 1);

    const auto inspection = SamFileUtility::inspect(path("sort_missing_eof.bam"));

    EXPECT_EQ(inspection.status, SamFileUtility::SamFileStatus::HasRecords);
    EXPECT_FALSE(inspection.hasMissingEof());
    EXPECT_EQ(inspection.recordCount, 1U);
}

TEST_F(SamFileUtilityTest, SortTruncatedInputFailsWithoutCreatingFinalOutput) {
    writeRecordBam(path("truncated.bam"));
    truncateMidBlock(path("truncated.bam"));

    ASSERT_EXIT(SamFileUtility::sortByQueryName(path("truncated.bam"), path("sorted.bam"), 1),
                testing::ExitedWithCode(EXIT_FAILURE), "truncated.bam");
    EXPECT_FALSE(fs::exists(path("sorted.bam")));
}

TEST_F(SamFileUtilityTest, MergeTruncatedInputFailsWithInputPath) {
    writeRecordBam(path("truncated.bam"));
    truncateMidBlock(path("truncated.bam"));

    ASSERT_EXIT(helper::mergeSamFiles({path("truncated.bam")}, path("merged.bam"),
                                      testReference()),
                testing::ExitedWithCode(EXIT_FAILURE), "truncated.bam");
}

TEST_F(SamFileUtilityTest, MergeWithoutReferenceOrUsableInputFailsClearly) {
    ASSERT_EXIT(helper::mergeSamFiles({}, path("merged.bam"), std::nullopt),
                testing::ExitedWithCode(EXIT_FAILURE),
                "without records or a SAM reference");
}

// NOLINTEND(readability-magic-numbers)
