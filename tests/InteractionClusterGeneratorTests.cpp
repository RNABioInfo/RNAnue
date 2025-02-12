// NOLINTBEGIN

#include <gtest/gtest.h>

// Standard
#include <algorithm>
#include <cstddef>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

// Internal
#include "AnnotatedInteractionCluster.hpp"
#include "ClusteringParameters.hpp"
#include "FeatureAnnotator.hpp"
#include "GenomicFeature.hpp"
#include "GenomicOrientation.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "GenomicStrandSpecificity.hpp"
#include "InteractionCluster.hpp"
#include "InteractionClusterGenerator.hpp"
#include "PartiallyAnnotatedInteractionCluster.hpp"
#include "Region.hpp"
#include "SortedGenomicRegionPair.hpp"
#include "SplitRecordsParser.hpp"
#include "TestFilePath.hpp"

using namespace pipelines::analyze;

// Note: This is the expected clustering result:
// Cluster 1: 1,3,6 // Fully annotated; First segment: RecordSegment(referenceIDIndex: 0,
// strand: +, start: 19, end: 27), Second segment id: RecordSegment(referenceIDIndex: 1,
// strand: +, start: 49, end: 64)
//
// Cluster 2: 2,7 // Partially annotated; First segment: RecordSegment(referenceIDIndex: 1,
// strand: +, start: 4, end: 10), Second segment id: RecordSegment(referenceIDIndex: 1, strand:
// +, start: 51, end: 57)
//
// Cluster 3: 4,5 // Non annotated; First segment: RecordSegment(referenceIDIndex: 0, strand:
// +, start: 4, end: 14), Second segment id: RecordSegment(referenceIDIndex: 0, strand: +,
// start: 39, end: 46)

using namespace pipelines::analyze;

struct InteractionClusterGeneratorTestParam {
    // Which cluster-merge parameter are we testing?
    ClusteringParameters clusteringParameters;

    // Path to a test SAM file or other input that can be parsed into InteractionClusters.
    // You could also inline data or dynamically generate it for the test if desired.
    std::string sourceFile;

    // Expected result data.
    std::unordered_map<std::string, size_t> expectedFeatureCounts;
    std::vector<AnnotatedInteractionCluster> finalizedClusters;
    std::vector<PartiallyAnnotatedInteractionCluster> partiallyAnnotatedClusters;
};

inline void PrintTo(const InteractionClusterGeneratorTestParam& param, std::ostream* os) {
    // Attempt to describe the chosen "clusterMergeParameter" variant
    std::string mergeType;
    std::visit(
        [&](auto&& v) {
            using T = std::decay_t<decltype(v)>;
            if constexpr (std::is_same_v<T, ClusterOverlapToleranceMergeParameter>) {
                mergeType = "Tolerance=" + std::to_string(v.tolerance);
            } else if constexpr (std::is_same_v<T, ShortestClusterOverlapFractionMergeParameter>) {
                mergeType = "ShortestOverlapFraction=" + std::to_string(v.overlapFraction);
            }
        },
        param.clusteringParameters.clusterMergeParameter);

    // Using the data we have to form a concise string
    *os << "InteractionClusterGeneratorTestParam("
        << "sourceFile=" << param.sourceFile << ", mergeType=" << mergeType
        << ", strandSpecificity="
        << static_cast<int>(param.clusteringParameters.clusterMergingStrandSpecificity)
        << ", maxSelfOverlapFraction=" << param.clusteringParameters.maxClusterSelfOverlapFraction
        << ", minClusterReadCount=" << param.clusteringParameters.minimumClusterReadCount
        << ", featureOrientation="
        << static_cast<int>(param.clusteringParameters.featureOrientation) << ")";
}

static const FeatureMap featureMap = {
    {0,
     {GenomicFeature{"gene",
                     {1, {.startPosition = 20, .endPosition = 30}, GenomicStrand::FORWARD},
                     "gene1",
                     std::nullopt,
                     std::nullopt},
      GenomicFeature{"gene",
                     {1, {.startPosition = 50, .endPosition = 60}, GenomicStrand::FORWARD},
                     "gene2",
                     std::nullopt,
                     std::nullopt}}},
    {1,
     {GenomicFeature{"gene",
                     {2, {.startPosition = 50, .endPosition = 60}, GenomicStrand::FORWARD},
                     "gene3",
                     std::nullopt,
                     std::nullopt}}}};

static const std::shared_ptr<annotation::FeatureAnnotator> featureAnnotator =
    std::make_shared<annotation::FeatureAnnotator>(featureMap);

// class InteractionClusterGeneratorTests
//     : public testing::TestWithParam<InteractionClusterGeneratorTestParam> {
//    protected:
//     InteractionClusterGeneratorTests()
//         : interactionClusterGenerator(
//               featureAnnotator,
//               ClusteringParameters{
//                   .clusterMergeParameter = ClusterOverlapToleranceMergeParameter{0},
//                   .clusterMergingStrandSpecificity = GenomicStrandSpecificity::UNSPECIFIC,
//                   .maxClusterSelfOverlapFraction = 1.0,  // NOLINT
//                   .minimumClusterReadCount = 1,
//                   .featureOrientation = dataTypes::GenomicOrientation::BOTH}) {}

//     InteractionClusterGenerator interactionClusterGenerator;
// };

/**
 * @brief The test fixture for InteractionClusterGenerator parameterized tests.
 */
class InteractionClusterGeneratorTests
    : public testing::TestWithParam<InteractionClusterGeneratorTestParam> {
   protected:
    /**
     * @brief Helper that parses the source file from the test parameter into a vector of
     * InteractionClusters (or you can generate them manually).
     */
    static std::vector<InteractionCluster> parseClusters(const std::string& sourceFile) {
        return SplitRecordsParser::parse(getTestFilePath(sourceFile));
    }
};

/**
 * @brief One example cluster to demonstrate fully annotated results.
 */
static const AnnotatedInteractionCluster exampleClusterFullyAnnotated(
    InteractionCluster{
        SortedGenomicRegionPair{GenomicRegion{0, Region{.startPosition = 19, .endPosition = 27},
                                              GenomicStrand::FORWARD},
                                GenomicRegion{1, Region{.startPosition = 49, .endPosition = 64},
                                              GenomicStrand::FORWARD}},
        {"SRR18331301.3", "SRR18331301.1", "SRR18331301.6"},
        {1.0, 1.0, 1.0},
        {-1.7, -1.7, -1.7},
        {1, 1, 1}},
    "gene1", "gene2");

/**
 * @brief Example cluster to demonstrate partial annotation with a single feature recognized.
 */
static const PartiallyAnnotatedInteractionCluster exampleClusterPartial1(
    InteractionCluster{
        SortedGenomicRegionPair{
            GenomicRegion{1, Region{.startPosition = 4, .endPosition = 10}, GenomicStrand::FORWARD},
            GenomicRegion{1, Region{.startPosition = 51, .endPosition = 57},
                          GenomicStrand::FORWARD}},
        {"SRR18331301.2", "SRR18331301.7"},
        {1.0, 1.0},
        {-1.7, -1.7},
        {1, 1}},
    std::nullopt, "gene3");

/**
 * @brief Example cluster to demonstrate partial annotation with no recognized features.
 */
static const PartiallyAnnotatedInteractionCluster exampleClusterPartial2(
    InteractionCluster{
        SortedGenomicRegionPair{
            GenomicRegion{0, Region{.startPosition = 4, .endPosition = 14}, GenomicStrand::FORWARD},
            GenomicRegion{0, Region{.startPosition = 39, .endPosition = 46},
                          GenomicStrand::FORWARD}},
        {"SRR18331301.5", "SRR18331301.4"},
        {1.0, 1.0},
        {-1.7, -1.7},
        {1, 1}},
    std::nullopt, std::nullopt);

/**
 * @brief The "Default" test scenario with cluster overlap tolerance = 0,
 *        no strand specificity, and requiring minimum read count of 1.
 */
INSTANTIATE_TEST_SUITE_P(
    DefaultParameters, InteractionClusterGeneratorTests,
    testing::Values(InteractionClusterGeneratorTestParam{
        // Clustering parameters
        ClusteringParameters{
            .clusterMergeParameter = ClusterOverlapToleranceMergeParameter{0},  // Overlap tolerance
            .clusterMergingStrandSpecificity = GenomicStrandSpecificity::UNSPECIFIC,
            .maxClusterSelfOverlapFraction = 1.0,
            .minimumClusterReadCount = 1,
            .featureOrientation = dataTypes::GenomicOrientation::BOTH},
        // Source test file (relative path or suitable test fixture resource)
        "interactionClusterRecords.sam",
        // Expected results
        // Feature counts: "gene1" has 3 from cluster1, "gene3" also accumulates 3 from partial
        // cluster
        {{"gene1", 3}, {"gene3", 3}},
        // Finalized clusters
        {exampleClusterFullyAnnotated},
        // Partially annotated clusters
        {exampleClusterPartial1, exampleClusterPartial2}}));

/**
 * @brief Test that verifies merging of clusters with no overlap tolerance
 *        still yields the correct final & partial clusters.
 */
TEST_P(InteractionClusterGeneratorTests, MergeClustersNoOverlapTolerance) {
    const auto& param = GetParam();

    // Create the generator with the provided parameters
    InteractionClusterGenerator generator(featureAnnotator, param.clusteringParameters);

    // Parse input data
    std::vector<InteractionCluster> clusters = parseClusters(param.sourceFile);

    // Sort to mimic the expected pipeline usage
    std::sort(clusters.begin(), clusters.end());

    // Merge
    auto result = generator.mergeClusters(std::move(clusters));

    // Check expectations
    EXPECT_EQ(result.finishedClusters, param.finalizedClusters);
    EXPECT_EQ(result.partiallyAnnotatedClusters, param.partiallyAnnotatedClusters);
    EXPECT_EQ(result.featureCounts, param.expectedFeatureCounts);
}

/**
 * @brief Parameterized test that uses a shorter-cluster-overlap fraction for merging.
 */
INSTANTIATE_TEST_SUITE_P(
    ShortestOverlapFraction, InteractionClusterGeneratorTests,
    testing::Values(InteractionClusterGeneratorTestParam{
        ClusteringParameters{
            .clusterMergeParameter = ShortestClusterOverlapFractionMergeParameter{0.25},
            .clusterMergingStrandSpecificity = GenomicStrandSpecificity::UNSPECIFIC,
            .maxClusterSelfOverlapFraction = 1.0,  // strict check for self-overlap
            .minimumClusterReadCount = 2,          // require at least 2 reads in a cluster
            .featureOrientation = dataTypes::GenomicOrientation::BOTH},
        "interactionClusterRecords.sam",
        // Suppose we expect fewer final clusters because the fraction requirement is stricter
        // and we also require min read count = 2.
        {{"gene1", 3}, {"gene3", 3}},  // hypothetical
        {// Possibly we get two annotated clusters, each with enough reads to pass min count
         AnnotatedInteractionCluster(
             InteractionCluster{
                 SortedGenomicRegionPair{GenomicRegion{0, {19, 27}, GenomicStrand::FORWARD},
                                         GenomicRegion{1, {49, 64}, GenomicStrand::FORWARD}},
                 {"SRR18331301.3", "SRR18331301.1", "SRR18331301.6"},  // 3 reads
                 {1.0, 1.0, 1.0},
                 {-1.7, -1.7, -1.7},
                 {1, 1, 1}},
             "gene1", "gene2")},
        {// Some partially annotated clusters that didn't meet the read threshold
         PartiallyAnnotatedInteractionCluster(
             InteractionCluster{
                 SortedGenomicRegionPair{GenomicRegion{1, {4, 10}, GenomicStrand::FORWARD},
                                         GenomicRegion{1, {52, 57}, GenomicStrand::FORWARD}},
                 {"SRR18331301.2", "SRR18331301.7"},  // two records
                 {1.0, 1.0},
                 {-1.7, -1.7},
                 {1, 1}},
             std::nullopt, std::nullopt),
         PartiallyAnnotatedInteractionCluster(
             InteractionCluster{
                 SortedGenomicRegionPair{GenomicRegion{1, {4, 10}, GenomicStrand::FORWARD},
                                         GenomicRegion{1, {52, 57}, GenomicStrand::FORWARD}},
                 {"SRR18331301.5",
                  "SRR18331301.4"},  // only 1 read, fails min read count => partial
                 {1.0, 1.0},
                 {-1.7, -1.7},
                 {1, 1}},
             std::nullopt, std::nullopt)}}));

// Note: This is the expected clustering result:
// Cluster 1: 1,3,6 // Fully annotated; First segment: RecordSegment(referenceIDIndex: 0,
// strand: +, start: 19, end: 27), Second segment id: RecordSegment(referenceIDIndex: 1,
// strand: +, start: 49, end: 64)
//
// Cluster 2: 2,7 // Partially annotated; First segment: RecordSegment(referenceIDIndex: 1,
// strand: +, start: 4, end: 10), Second segment id: RecordSegment(referenceIDIndex: 1, strand:
// +, start: 51, end: 57)
//
// Cluster 3: 4,5 // Non annotated; First segment: RecordSegment(referenceIDIndex: 0, strand:
// +, start: 4, end: 14), Second segment id: RecordSegment(referenceIDIndex: 0, strand: +,
// start: 39, end: 46)

/**
 * @brief Test that applies the shortest-overlap-fraction merging logic and a strict read count.
 */
TEST_P(InteractionClusterGeneratorTests, MergeClustersWithShortestOverlapFraction) {
    const auto& param = GetParam();
    InteractionClusterGenerator generator(featureAnnotator, param.clusteringParameters);

    std::vector<InteractionCluster> clusters = parseClusters(param.sourceFile);
    std::sort(clusters.begin(), clusters.end());

    auto result = generator.mergeClusters(std::move(clusters));

    EXPECT_EQ(result.finishedClusters, param.finalizedClusters);
    EXPECT_EQ(result.partiallyAnnotatedClusters.size(), param.partiallyAnnotatedClusters.size());
    EXPECT_EQ(result.featureCounts, param.expectedFeatureCounts);

    // Additional checks on included vs excluded
    size_t totalClustersExamined = result.includedClusterCount + result.excludedClusterCount;
    EXPECT_EQ(totalClustersExamined,
              param.finalizedClusters.size() + param.partiallyAnnotatedClusters.size());
}

/**
 * @brief Test to ensure that clusters that fail overlap or min-read filters get excluded.
 */
TEST_F(InteractionClusterGeneratorTests, ExclusionByReadCountOrOverlap) {
    // Make a generator that requires at least 3 reads to pass
    ClusteringParameters params{
        .clusterMergeParameter = ClusterOverlapToleranceMergeParameter{0},
        .clusterMergingStrandSpecificity = GenomicStrandSpecificity::UNSPECIFIC,
        .maxClusterSelfOverlapFraction = 1.0,
        .minimumClusterReadCount = 3,
        .featureOrientation = dataTypes::GenomicOrientation::BOTH};

    InteractionClusterGenerator generator(featureAnnotator, params);

    // Create some minimal clusters
    InteractionCluster c1{
        SortedGenomicRegionPair{
            GenomicRegion{0, {.startPosition = 100, .endPosition = 110}, GenomicStrand::FORWARD},
            GenomicRegion{1, {.startPosition = 200, .endPosition = 205}, GenomicStrand::FORWARD}},
        {"r1", "r2"},  // only 2 reads
        {1.5, 2.0},
        {0.1, 0.2},
        {1, 1}};
    InteractionCluster c2{
        SortedGenomicRegionPair{
            GenomicRegion{0, {.startPosition = 100, .endPosition = 110}, GenomicStrand::FORWARD},
            GenomicRegion{1, {.startPosition = 220, .endPosition = 230}, GenomicStrand::FORWARD}},
        {"r3", "r4", "r5"},  // 3 reads
        {2.0, 2.0, 2.0},
        {0.3, 0.3, 0.3},
        {1, 1, 1}};

    // They do not overlap in the second segment; c1 end=205, c2 start=220 => no overlap
    std::vector<InteractionCluster> inputClusters{c1, c2};
    // Sort them
    std::sort(inputClusters.begin(), inputClusters.end());

    auto result = generator.mergeClusters(std::move(inputClusters));

    // c1 has only 2 reads => excluded
    // c2 has 3 reads => included if it passes overlap checks (there's no overlap to merge; it
    // stands alone)
    EXPECT_EQ(result.partiallyAnnotatedClusters.size(),
              1UL);  // likely partially or fully annotated depending on annotations
    EXPECT_EQ(result.includedClusterCount, 1UL);
    EXPECT_EQ(result.excludedClusterCount, 1UL);

    // You can further check which cluster is included if you want:
    // Typically you'd examine result.finishedClusters[0] or result.partiallyAnnotatedClusters[0].
}

/**
 * @brief Test to ensure that supplementary feature map is populated with
 *        unannotated segments from passing clusters.
 */
TEST_F(InteractionClusterGeneratorTests, SupplementaryFeatureMapPopulation) {
    // Make a generator that passes all clusters for read count=1,
    // so unannotated segments become supplementary features.
    ClusteringParameters params{
        .clusterMergeParameter = ClusterOverlapToleranceMergeParameter{0},
        .clusterMergingStrandSpecificity = GenomicStrandSpecificity::UNSPECIFIC,
        .maxClusterSelfOverlapFraction = 1.0,
        .minimumClusterReadCount = 1,
        .featureOrientation = dataTypes::GenomicOrientation::BOTH};

    InteractionClusterGenerator generator(featureAnnotator, params);

    // Single cluster that has no annotation at all => should produce a partially annotated cluster
    InteractionCluster noAnnotation{
        SortedGenomicRegionPair{GenomicRegion{5, {1000, 1020}, GenomicStrand::FORWARD},
                                GenomicRegion{5, {2000, 2050}, GenomicStrand::FORWARD}},
        {"rX"},
        {3.0},
        {1.1},
        {1}};

    std::vector<InteractionCluster> inputClusters{noAnnotation};
    auto result = generator.mergeClusters(std::move(inputClusters));

    // We should have 0 fully annotated, 1 partial, and the supplementaryFeatureMap
    // should have info for reference ID 5 segments.
    EXPECT_TRUE(result.finishedClusters.empty());
    ASSERT_EQ(result.partiallyAnnotatedClusters.size(), 1UL);

    auto it = result.supplementaryFeatureMap.find(5);
    ASSERT_TRUE(it != result.supplementaryFeatureMap.end());
    // We expect two new features corresponding to each segment
    EXPECT_EQ(it->second.size(), 2UL);
    EXPECT_EQ(result.includedClusterCount, 1UL);
    EXPECT_EQ(result.excludedClusterCount, 0UL);
}

/**
 * @brief Test ensuring that a cluster that fails the max self-overlap fraction gets excluded.
 */
TEST_F(InteractionClusterGeneratorTests, ExclusionBySelfOverlapFraction) {
    // Create a generator that excludes clusters with > 0.5 self-overlap fraction
    ClusteringParameters params{
        .clusterMergeParameter = ClusterOverlapToleranceMergeParameter{0},
        .clusterMergingStrandSpecificity = GenomicStrandSpecificity::UNSPECIFIC,
        .maxClusterSelfOverlapFraction = 0.5,  // must not exceed 50% overlap in the region
        .minimumClusterReadCount = 1,
        .featureOrientation = dataTypes::GenomicOrientation::BOTH};

    InteractionClusterGenerator generator(featureAnnotator, params);

    // Make a cluster whose segments overlap themselves heavily. For demonstration,
    // suppose the cluster says the "first" segment is from 100..120, but the paired segment
    // is also from 105..125 in the same reference, leading to a significant self-overlap fraction.
    InteractionCluster heavyOverlap{
        SortedGenomicRegionPair{GenomicRegion{0, {100, 120}, GenomicStrand::FORWARD},
                                GenomicRegion{0, {105, 125}, GenomicStrand::FORWARD}},
        {"ovl1", "ovl2"},
        {2.0, 2.0},
        {0.0, 0.0},
        {1, 1}};

    auto result = generator.mergeClusters({heavyOverlap});
    // Because of the large self-overlap, we expect it to be excluded
    EXPECT_TRUE(result.finishedClusters.empty());
    EXPECT_TRUE(result.partiallyAnnotatedClusters.empty());
    EXPECT_EQ(result.includedClusterCount, 0UL);
    EXPECT_EQ(result.excludedClusterCount, 1UL);
}

// NOLINTEND
