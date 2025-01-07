#include <gtest/gtest.h>

#include <seqan3/alignment/cigar_conversion/alignment_from_cigar.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/utility/type_list/type_list.hpp>
#include <sstream>

#include "FeatureAnnotator.hpp"
#include "ParseSamRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"
#include "SplitRecordsSplicingEvaluator.hpp"

using namespace seqan3::literals;

struct IsSplicedTestParam {
    SplitRecords splitRecords;
    bool isSpliced;
};

class EvaluatedSplitRecordsTests : public testing::TestWithParam<IsSplicedTestParam> {
   protected:
    EvaluatedSplitRecordsTests() : featureAnnotator(featureMap) {};

    const std::deque<std::string> referenceIDs = {"chromosome1"};

    const dataTypes::FeatureMap featureMap = {{"chromosome1",
                                               {{.referenceID = "chromosome1",
                                                 .type = "exon",
                                                 .startPosition = 1,
                                                 .endPosition = 10,
                                                 .strand = dataTypes::GenomicStrand::FORWARD,
                                                 .id = "exon1",
                                                 .groupID = "gene1",
                                                 .geneName = std::nullopt},
                                                {.referenceID = "chromosome1",
                                                 .type = "exon",
                                                 .startPosition = 20,
                                                 .endPosition = 30,
                                                 .strand = dataTypes::GenomicStrand::FORWARD,
                                                 .id = "exon2",
                                                 .groupID = "gene1",
                                                 .geneName = std::nullopt},
                                                {.referenceID = "chromosome1",
                                                 .type = "exon",
                                                 .startPosition = 40,
                                                 .endPosition = 50,
                                                 .strand = dataTypes::GenomicStrand::FORWARD,
                                                 .id = "exon3",
                                                 .groupID = "gene1",
                                                 .geneName = std::nullopt}}}};

    annotation::FeatureAnnotator featureAnnotator;
};

const auto noSpliceRaw = R"(
@HD     VN:1.6
@SQ     SN:chromosome1 LN:100
SRR18331301.231	0	chromosome1	5	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.232	0	chromosome1	10	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

const auto spliceRaw = R"(
@HD     VN:1.6
@SQ     SN:chromosome1 LN:100
SRR18331301.231	0	chromosome1	5	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.232	0	chromosome1	20	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

const auto noSpliceInBetweenExonRaw = R"(
@HD     VN:1.6
@SQ     SN:chromosome1 LN:100
SRR18331301.231	0	chromosome1	5	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.232	0	chromosome1	41	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

void PrintTo(const IsSplicedTestParam& param, std::ostream* os) {
    *os << "IsSplicedTestParam{.isSpliced = " << param.isSpliced << "}";
};

TEST_P(EvaluatedSplitRecordsTests, IsSplicedSplitRecord) {
    const IsSplicedTestParam& param = GetParam();

    const SplitRecordsEvaluationParameters::SplicingParameters splicingParameters = {
        .baseParameters =
            SplitRecordsEvaluationParameters::BaseParameters{
                .minComplementarity = 0.9,
                .minComplementarityFraction = 0.9,
                .mfeThreshold = 10,
                .includeWobbleBasePairsInCrosslinkingSites = true},
        .orientation = annotation::Orientation::BOTH,
        .splicingTolerance = 0,
        .featureAnnotator = &featureAnnotator};

    const auto isSpliced = SplitRecordsSplicingEvaluator::isSplicedSplitRecord(
        param.splitRecords, referenceIDs, splicingParameters);

    EXPECT_EQ(isSpliced, param.isSpliced);
};

INSTANTIATE_TEST_SUITE_P(
    Default, EvaluatedSplitRecordsTests,
    testing::Values(
        IsSplicedTestParam{.splitRecords = parseSamRecords(noSpliceRaw), .isSpliced = false},
        IsSplicedTestParam{.splitRecords = parseSamRecords(spliceRaw), .isSpliced = true},
        IsSplicedTestParam{.splitRecords = parseSamRecords(noSpliceInBetweenExonRaw),
                           .isSpliced = false}));
