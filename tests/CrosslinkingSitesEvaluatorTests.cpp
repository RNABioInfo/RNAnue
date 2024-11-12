#include <gtest/gtest.h>

#include "gtest/gtest.h"

// Standard
#include <vector>

// seqan3
#include <seqan3/alphabet/structure/dot_bracket3.hpp>

// Internal
#include "CrosslinkingSitesEvaluator.hpp"

using namespace pipelines::detect;

struct CrosslinkingSitesEvaluatorTestParams {
    seqan3::dna5_vector sequence1;
    seqan3::dna5_vector sequence2;
    std::vector<seqan3::dot_bracket3> dotbracket;

    CrosslinkingSitesEvaluator::Result expected;
};

class CrosslinkingSitesEvaluatorTest
    : public ::testing::TestWithParam<CrosslinkingSitesEvaluatorTestParams> {};

TEST_P(CrosslinkingSitesEvaluatorTest, Base) {
    CrosslinkingSitesEvaluatorTestParams param = GetParam();
    const auto result =
        CrosslinkingSitesEvaluator::evaluate(param.sequence1, param.sequence2, param.dotbracket);
    EXPECT_EQ(result, param.expected);
}

const CrosslinkingSitesEvaluatorTestParams params1{
    .sequence1 = "CCCUACCC"_dna5,
    .sequence2 = "GGGUAGGG"_dna5,
    .dotbracket = {"((((((((&))))))))"_db3},
    .expected = CrosslinkingSitesEvaluator::Result{{{}, {}}, {{3, 3}}, "((((((((&))))))))"}};

const CrosslinkingSitesEvaluatorTestParams params2{
    .sequence1 = "UUAACAAACCACCUGCAU"_dna5,
    .sequence2 = "GUAUCCCGUAGGGGAGUC"_dna5,
    .dotbracket = {"........((.(((((..&.......)))))))...."_db3},
    .expected = CrosslinkingSitesEvaluator::Result{
        {{}, {}}, {{13, 8}}, "........((.(((((..&.......)))))))...."}};

const CrosslinkingSitesEvaluatorTestParams params3{
    .sequence1 = "ACCCGACAAGGAAUUUCGC"_dna5,
    .sequence2 = "UGCGCCCAUUGUGCAAU"_dna5,
    .dotbracket = {"..((.....))........&(((((.....))))).."_db3},
    .expected = CrosslinkingSitesEvaluator::Result{
        {{}, {{2, 13}}}, {}, "..((.....))........&(((((.....))))).."}};

const CrosslinkingSitesEvaluatorTestParams params4{
    .sequence1 = "CAGAGCCGCUGCUUUGA"_dna5,
    .sequence2 = "AGCCAAUCCGCUAGACGCU"_dna5,
    .dotbracket = {"((((((....)))))).&(((......)))......."_db3},
    .expected = CrosslinkingSitesEvaluator::Result{
        {{{5, 11}}, {{2, 10}}}, {}, "((((((....)))))).&(((......)))......."}};

const CrosslinkingSitesEvaluatorTestParams params5{
    .sequence1 = "UCUGAAGCAGACUGCAAAG"_dna5,
    .sequence2 = "GUUAGAGCGUACGCCUG"_dna5,
    .dotbracket = {"(((((.(((...)))....&.)))))((....))..."_db3},
    .expected = CrosslinkingSitesEvaluator::Result{
        {{{7, 14}}, {{7, 13}}}, {{2, 2}}, "(((((.(((...)))....&.)))))((....))..."}};

INSTANTIATE_TEST_SUITE_P(Default, CrosslinkingSitesEvaluatorTest,
                         testing::Values(params1, params2, params3, params4, params5));
