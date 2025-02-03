# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------
include(FetchContent)
FetchContent_Declare(
    googletest
    DOWNLOAD_EXTRACT_TIMESTAMP true
    URL https://github.com/google/googletest/archive/1d17ea141d2c11b8917d2c7d029f1c4e2b9769b2.zip
)

FetchContent_MakeAvailable(googletest)

enable_testing ()

add_executable (tests   ${CMAKE_CURRENT_SOURCE_DIR}/InteractionClusterTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/InteractionClusterGeneratorTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/ParallelInteractionClusterGeneratorTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/FeatureAnnotatorTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/DetectTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/SplitRecordsSplicingEvaluatorTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/SplitRecordsTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/AsyncSplitRecordGroupBufferTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/StatisticEvaluatorTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/DeduplicatorTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/CrosslinkingSitesEvaluatorTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/IITreeTests.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/GenomicFeatureTreeMergerTests.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/FeatureAnnotator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/FeatureParser.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/FeatureWriter.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/GenomicFeatureTreeMerger.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/utility/SamFileUtility.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/utility/SequenceFileUtility.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/dataTypes/GenomicRegion.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/dataTypes/RecordFragment.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/dataTypes/GenomicOrientation.cpp


                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/preprocess/Deduplicator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/preprocess/DeduplicatorBySequenceSingleEnd.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/preprocess/DeduplicatorBySequencePairedEnd.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/CrosslinkingSitesEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/Detect.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/DetectData.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/SplitRecordsComplementarityEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/SplitRecordsEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/SplitRecordsHybridizationEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/SplitRecordsSplicingEvaluator.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/Analyze.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/AnalyzeData.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/InteractionCluster.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/InteractionClusterGenerator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/ParallelInteractionClusterGenerator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/InteractionsWriter.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/StatisticEvaluator.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/Utility.cpp
)

# # # htslib
target_link_libraries (tests PRIVATE ${HTSlib_LIBRARIES})

# # # ViennaRNA
target_link_libraries (tests PRIVATE ${VIENNA_RNA_LIBRARY})

target_include_directories (tests PRIVATE   ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                            ${CMAKE_CURRENT_SOURCE_DIR}/../include/annotation
                                            ${CMAKE_CURRENT_SOURCE_DIR}/../include/pipelines/detect
                                            ${CMAKE_CURRENT_SOURCE_DIR}/../include/pipelines/dataTypes
                                            ${CMAKE_CURRENT_SOURCE_DIR}/../include/utility
                                            ${VIENNA_RNA_INCLUDE_DIR}
                                            ${VIENNA_RNA_INCLUDE_DIR}/ViennaRNA
                                            ${HTSlib_INCLUDE_DIRS}
                                            ${HTSlib_INCLUDE_DIRS}/htslib
                                            ${Boost_INCLUDE_DIR})

target_link_libraries(tests     PRIVATE
                                GTest::gtest_main
                                seqan3::seqan3
                                segemehl
                                ${Boost_LIBRARIES})


include(GoogleTest)
gtest_discover_tests(tests)
