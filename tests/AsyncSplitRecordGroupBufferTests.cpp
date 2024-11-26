#include <gtest/gtest.h>

// Standard
#include <cstddef>
#include <filesystem>
#include <future>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>

// seqan3
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>

// Internal
#include "AsyncSplitRecordGroupBuffer.hpp"
#include "TestFilePath.hpp"

using namespace seqan3::literals;

TEST(AsyncSplitRecordGroupBufferTest, SingleThreaded) {
    seqan3::sam_file_input fin{getTestFilePath("splitRecords.bam"), dataTypes::SamFieldIDs{}};

    auto v = fin | AsyncSplitRecordGroupBuffer(2);

    size_t groupCount = 0;

    for (auto& group : v) {
        ASSERT_NE(group.size(), 0ul);

        std::set<std::string> doubleIDs = {"SRR18331301.2", "SRR18331301.6", "SRR18331301.7"};
        std::set<std::string> tripleIDs = {"SRR18331301.16", "SRR18331301.26"};

        if (group.front().id() == "SRR18331301.1") {
            EXPECT_EQ(group.size(), 1ul);
        } else if (doubleIDs.contains(group.front().id())) {
            EXPECT_EQ(group.size(), 2ul);
        } else if (tripleIDs.contains(group.front().id())) {
            EXPECT_EQ(group.size(), 3ul);
        } else {
            throw std::logic_error("File should not contain groups of size: " +
                                   std::to_string(group.size()));
        }

        groupCount++;
    }

    EXPECT_EQ(groupCount, 6ul);
};

TEST(AsyncSplitRecordGroupBufferTest, Multithreaded) {
    seqan3::sam_file_input fin{getTestFilePath("splitRecords.bam"), dataTypes::SamFieldIDs{}};

    auto asyncInputBuffer = fin | AsyncSplitRecordGroupBuffer(2);

    auto worker = [&asyncInputBuffer]() -> size_t {
        size_t count = 0;

        for (auto& group : asyncInputBuffer) {
            std::cout << "Thread ID: " << std::this_thread::get_id() << "\n";
            std::cout << "GROUP ID: " << group.front().id() << std::endl;
            count++;

            std::set<std::string> doubleIDs = {"SRR18331301.2", "SRR18331301.6", "SRR18331301.7"};
            std::set<std::string> tripleIDs = {"SRR18331301.16", "SRR18331301.26"};

            if (group.front().id() == "SRR18331301.1") {
                EXPECT_EQ(group.size(), 1ul);
            } else if (doubleIDs.contains(group.front().id())) {
                EXPECT_EQ(group.size(), 2ul);
            } else if (tripleIDs.contains(group.front().id())) {
                EXPECT_EQ(group.size(), 3ul);
            } else {
                throw std::logic_error("File should not contain groups of size: " +
                                       std::to_string(group.size()));
            }
        }

        return count;
    };

    auto f0 = std::async(std::launch::async, worker);
    auto f1 = std::async(std::launch::async, worker);

    EXPECT_EQ(f0.get() + f1.get(), 6ul);
}
