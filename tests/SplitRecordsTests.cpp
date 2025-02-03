#include <gtest/gtest.h>

// Standard
#include <algorithm>
#include <functional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

// Internal
#include "ParseSamRecords.hpp"
#include "SplitRecords.hpp"

using namespace dataTypes;

struct SplitRecordsTestParams {
    SplitRecordsTestParams(std::vector<SplitRecords> splitRecords,
                           std::vector<std::string> expectedBackRecordIDOrder)
        : splitRecords(std::move(splitRecords)),
          expectedBackRecordIDOrder(std::move(expectedBackRecordIDOrder)) {}
    const std::vector<SplitRecords> splitRecords;
    const std::vector<std::string> expectedBackRecordIDOrder;
};

class SplitRecordsTests : public ::testing::TestWithParam<SplitRecordsTestParams> {};

TEST_P(SplitRecordsTests, IsSortedFromBackToFront) {
    const SplitRecordsTestParams& param = GetParam();
    std::vector<SplitRecords> splitRecordGroups = param.splitRecords;
    const std::vector<std::string>& expectedBackRecordIDOrder = param.expectedBackRecordIDOrder;

    EXPECT_EQ(expectedBackRecordIDOrder.size(), splitRecordGroups.size());

    std::ranges::sort(splitRecordGroups, std::greater());

    std::vector<std::string> backRecordIDOrder;
    backRecordIDOrder.reserve(splitRecordGroups.size());
    for (const auto& splitRecords : splitRecordGroups) {
        backRecordIDOrder.push_back(splitRecords.back().id());
    }

    EXPECT_EQ(expectedBackRecordIDOrder, backRecordIDOrder);
}

const auto splitRecords1 = R"(@HD	VN:1.6
@SQ	SN:chromosome1	LN:100
@SQ	SN:chromosome2	LN:100
SRR18331301.231	0	chromosome1	0	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.232	0	chromosome1	40	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

const auto splitRecords2 = R"(@HD	VN:1.6
@SQ	SN:chromosome1	LN:100
@SQ	SN:chromosome2	LN:100
SRR18331301.233	0	chromosome1	0	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.234	0	chromosome1	50	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

const auto splitRecords3 = R"(@HD	VN:1.6
@SQ	SN:chromosome1	LN:100
@SQ	SN:chromosome2	LN:100
SRR18331301.235	0	chromosome1	0	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.236	0	chromosome1	50	20	7M	*	0	0	ATCGCGT	@@@@@@@	AS:i:0	XS:i:0
)";

const auto splitRecords4 = R"(@HD	VN:1.6
@SQ	SN:chromosome1	LN:100
@SQ	SN:chromosome2	LN:100
SRR18331301.230	0	chromosome2	40	20	7M	*	0	0	ATCGCGT	@@@@@@@	AS:i:0	XS:i:0
SRR18331301.229	0	chromosome1	0	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

void PrintTo(const SplitRecordsTestParams& param, std::ostream* outStream) {
    *outStream << "SplitRecordsTestParams{" << "\n" << "splitRecords back IDs: \n";

    for (const auto& splitRecords : param.splitRecords) {
        *outStream << "\t" << splitRecords.back().id()
                   << " reference id: " << splitRecords.back().reference_id().value_or(-1) << "\n";
    }

    *outStream << "expectedBackRecordIDOrder: \n";

    for (const auto& recordID : param.expectedBackRecordIDOrder) {
        *outStream << recordID << "\n";
    }

    *outStream << "}" << "\n";
};

const SplitRecordsTestParams sameChromosomeTestParams =
    SplitRecordsTestParams{{parseSamRecords(splitRecords2), parseSamRecords(splitRecords3),
                            parseSamRecords(splitRecords1)},
                           {"SRR18331301.236", "SRR18331301.234", "SRR18331301.232"}};

const SplitRecordsTestParams differentChromosomeTestParams = SplitRecordsTestParams{
    {parseSamRecords(splitRecords2), parseSamRecords(splitRecords4), parseSamRecords(splitRecords3),
     parseSamRecords(splitRecords1)},
    {"SRR18331301.230", "SRR18331301.236", "SRR18331301.234", "SRR18331301.232"}};

INSTANTIATE_TEST_SUITE_P(Default, SplitRecordsTests,
                         testing::Values(sameChromosomeTestParams, differentChromosomeTestParams));
