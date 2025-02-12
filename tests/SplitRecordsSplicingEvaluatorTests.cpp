#include <gtest/gtest.h>

// seqan3
#include <memory>
#include <optional>
#include <ostream>
#include <seqan3/alignment/cigar_conversion/alignment_from_cigar.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

// Internal
#include "FeatureAnnotator.hpp"
#include "GenomicFeature.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "ParseSamRecords.hpp"
#include "Region.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsEvaluationParameters.hpp"
#include "SplitRecordsSplicingEvaluator.hpp"

using dataTypes::FeatureMap;

using namespace seqan3::literals;

struct IsSplicedTestParam {
    SplitRecords splitRecords;
    bool isSpliced;
    bool allowAltsplice;
    dataTypes::GenomicOrientation annotationOrientatation;
};

class SplitRecordSplicingEvaluatorTests : public testing::TestWithParam<IsSplicedTestParam> {
   protected:
    SplitRecordSplicingEvaluatorTests()
        : featureAnnotator(std::make_shared<const FeatureAnnotator>(featureMap)) {};

    const FeatureMap featureMap = {
        {0,
         {GenomicFeature{"exon",
                         GenomicRegion{0, Region{.startPosition = 1, .endPosition = 9},
                                       GenomicStrand::FORWARD},
                         "exon1", "gene1", std::nullopt},
          {"exon",
           GenomicRegion{0, Region{.startPosition = 19, .endPosition = 30}, GenomicStrand::FORWARD},
           "exon2", "gene1", std::nullopt},
          {"exon",
           GenomicRegion{0, Region{.startPosition = 39, .endPosition = 50}, GenomicStrand::FORWARD},
           "exon3", "gene1", std::nullopt}}}};

    std::shared_ptr<const FeatureAnnotator> featureAnnotator;
};

const auto noSpliceRaw = R"(
@HD     VN:1.6
@SQ     SN:chromosome1 LN:100
SRR18331301.231	0	chromosome1	5	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.232	0	chromosome1	15	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

const auto spliceRaw = R"(
@HD     VN:1.6
@SQ     SN:chromosome1 LN:100
SRR18331301.231	0	chromosome1	5	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.232	0	chromosome1	20	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

const auto noSpliceRevStrand = R"(
@HD     VN:1.6
@SQ     SN:chromosome1 LN:100
SRR18331301.231	16	chromosome1	5	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.232	16	chromosome1	20	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

const auto noSpliceInBetweenExonRaw = R"(
@HD     VN:1.6
@SQ     SN:chromosome1 LN:100
SRR18331301.231	0	chromosome1	5	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
SRR18331301.232	0	chromosome1	40	20	5M	*	0	0	ATCGC	@@@@@	AS:i:0	XS:i:0
)";

void PrintTo(const IsSplicedTestParam& param, std::ostream* ostream) {
    *ostream << "IsSplicedTestParam{.isSpliced = " << param.isSpliced << "}";
};

TEST_P(SplitRecordSplicingEvaluatorTests, IsSplicedSplitRecord) {
    const IsSplicedTestParam& param = GetParam();

    const SplitRecordsEvaluationParameters::SplicingParameters splicingParameters = {
        .baseParameters =
            SplitRecordsEvaluationParameters::BaseParameters{
                .minComplementarity = 0.9,
                .minComplementarityFraction = 0.9,
                .mfeThreshold = 10,
                .includeWobbleBasePairsInCrosslinkingSites = true},
        .orientation = param.annotationOrientatation,
        .splicingTolerance = 0,
        .allowAlternativeSplicing = param.allowAltsplice,
        .featureAnnotator = featureAnnotator};

    const auto isSpliced =
        SplitRecordsSplicingEvaluator::isSplicedSplitRecord(param.splitRecords, splicingParameters);

    EXPECT_EQ(isSpliced, param.isSpliced);
};

INSTANTIATE_TEST_SUITE_P(
    Default, SplitRecordSplicingEvaluatorTests,
    testing::Values(IsSplicedTestParam{.splitRecords = parseSamRecords(noSpliceRaw),
                                       .isSpliced = false,
                                       .allowAltsplice = false,
                                       .annotationOrientatation = GenomicOrientation::SAME},
                    IsSplicedTestParam{.splitRecords = parseSamRecords(spliceRaw),
                                       .isSpliced = true,
                                       .allowAltsplice = false,
                                       .annotationOrientatation = GenomicOrientation::SAME},
                    IsSplicedTestParam{.splitRecords = parseSamRecords(noSpliceRevStrand),
                                       .isSpliced = false,
                                       .allowAltsplice = false,
                                       .annotationOrientatation = GenomicOrientation::SAME},
                    IsSplicedTestParam{.splitRecords = parseSamRecords(noSpliceRevStrand),
                                       .isSpliced = true,
                                       .allowAltsplice = false,
                                       .annotationOrientatation = GenomicOrientation::BOTH},
                    IsSplicedTestParam{.splitRecords = parseSamRecords(noSpliceRevStrand),
                                       .isSpliced = true,
                                       .allowAltsplice = false,
                                       .annotationOrientatation = GenomicOrientation::OPPOSITE},
                    IsSplicedTestParam{.splitRecords = parseSamRecords(noSpliceInBetweenExonRaw),
                                       .isSpliced = false,
                                       .allowAltsplice = false,
                                       .annotationOrientatation = GenomicOrientation::SAME},
                    IsSplicedTestParam{.splitRecords = parseSamRecords(noSpliceInBetweenExonRaw),
                                       .isSpliced = true,
                                       .allowAltsplice = true,
                                       .annotationOrientatation = GenomicOrientation::SAME}));
