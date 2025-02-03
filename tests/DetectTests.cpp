#include <gtest/gtest.h>

#include "Detect.hpp"

// class DetectTest : public testing::Test {
//    protected:
//     DetectTest() {
//         detect = Detect(po::variables_map{
//             "minlen" : 15,
//             "cmplmin" : 0.8,
//             "sitelenratio" : 0.5,
//             "mapqmin" : 10,
//             "exclclipping" : true,
//             "annotationorientation" : annotation::Orientation::BOTH,
//             "features" : "",
//             "featuretypes" : "gene"
//         });
//     }

//     Detect detect;
// };

// TEST_F(DetectTest, TestEmptyInput) {
//     std::vector<SamRecord> emptyRecords;
//     auto result = detect.constructSplitRecords(emptyRecords);
//     EXPECT_EQ(result, std::nullopt);
// }

// TEST_F(DetectTest, TestValidSingleRecord) {
//     using namespace seqan3::literals;

//     auto sam_file_raw = R"( @HD     VN:1.0  SO:queryname    SS:queryname:natural
//                             @SQ     SN:NC_002505.1  LN:2961149
//                             @SQ     SN:NC_002506.1  LN:1072315
//                             @RG     ID:A1   SM:sample1      LB:library1     PU:unit1 PL:illumina
//                             @PG     ID:segemehl     VN:0.3.4        CL:-S -A 90 -U 15 -W 80 -Z 15
//                             -t 8 -m 15 -i
//                             /Users/christopherphd/Documents/projects/RNAnue_dev/test_data/RILseq/results_develop_03/GCF_000006745.1_ASM674v1_genomic.idx
//                             -d
//                             /Users/christopherphd/Documents/projects/RNAnue_dev/test_data/vibrio_cholerae/ncbi_dataset/data/GCF_000006745.1/GCF_000006745.1_ASM674v1_genomic.fna
//                             -q
//                             /Users/christopherphd/Documents/projects/RNAnue_dev/test_data/RILseq/results_develop_03/preprocess/trtms/wt_expenential/SRR18331301_preproc.fastq
//                             -o
//                             /Users/christopherphd/Documents/projects/RNAnue_dev/test_data/RILseq/results_develop_03/align/trtms/wt_expenential/SRR18331301_preproc_matched.sam
//                             SRR18331301.232 0       NC_002505.1     2678474 29      36=23N39= *
//                             0       0
//                             TTTGCAAAGGGTTGGTAAGTCGGGATGACCCCCTAGTAGTATTCGTCCGAGGCGCTACCTAAATAGCTTTCGGGG
//                             GGGIIIIIIIIIGIGGGGIGGGGGGGIIGIIIIIIIIIIGIIIIIIIIIIIIIIIIIGIIIIIIIIIIIIGGGGA
//                             XS:A:+  YQ:A:P  HI:i:0  NH:i:1 NM:i:0   MD:Z:75 XM:B:I,0,0
//                             XL:B:I,36,39    XX:i:1  XY:i:75 XI:i:0  XH:i:2  XJ:i:2  RG:Z:A1
//                             YZ:Z:0)";

//     seqan3::sam_file_input fin{std::istringstream{sam_file_raw}, dtp::sam_field_ids{}};

//     SamRecord validRecord{fin.front()};
//     auto result = detect.constructSplitRecords({validRecord});
//     ASSERT_TRUE(result.has_value());
//     EXPECT_EQ(result->size(), 2);  // Expect at least one split record
// }

// TEST(DetectTest, TestInvalidCigarString) {
//     SamRecord invalidCigarRecord{/* Populate with data including an invalid cigar string */};
//     auto result = Detect::constructSplitRecords({invalidCigarRecord});
//     EXPECT_EQ(result, std::nullopt);
// }

// TEST(DetectTest, TestSoftClippingExclusion) {
//     SamRecord softClippedRecord{/* Populate with data including soft clipping */};
//     // Assuming excludeSoftClipping is a static variable or can be set in some way for the test
//     Detect::excludeSoftClipping = true;
//     auto result = Detect::constructSplitRecords({softClippedRecord});
//     // Validate based on expected behavior with soft clipping excluded
// }

// TEST(DetectTest, TestExpectedSplitRecordsMismatch) {
//     SamRecord mismatchRecord{/* Populate with data that should result in a mismatch between
//     expected and actual splits */};
//     auto result = Detect::constructSplitRecords({mismatchRecord});
//     EXPECT_EQ(result, std::nullopt);
// }

// TEST(DetectTest, TestMultipleRecords) {
//     std::vector<SamRecord> multipleRecords{
//         /* Populate with multiple records that should be split correctly */};
//     auto result = Detect::constructSplitRecords(multipleRecords);
//     ASSERT_TRUE(result.has_value());
//     // Additional checks based on expected outcome
// }
