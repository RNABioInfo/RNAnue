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

                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/FeatureAnnotator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/FeatureParser.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/FeatureWriter.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/dataTypes/GenomicRegion.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/preprocess/Deduplicator.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/CrosslinkingSitesEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/Detect.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/SplitRecordsComplementarityEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/SplitRecordsEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/SplitRecordsHybridizationEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/detect/SplitRecordsSplicingEvaluator.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/Analyze.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/RecordFragment.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/InteractionCluster.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/InteractionClusterGenerator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/ParallelInteractionClusterGenerator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/InteractionsWriter.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/StatisticEvaluator.cpp


                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/Utility.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/../src/CooptimalPairwiseAligner.cpp
)

# # # htslib
target_link_libraries (tests PRIVATE ${HTSlib_LIBRARIES})

# # # ViennaRNA
target_link_libraries (tests PRIVATE ${VIENNA_RNA_LIBRARY})

target_include_directories (tests PRIVATE   ${CMAKE_CURRENT_SOURCE_DIR}/../include
                                            ${CMAKE_CURRENT_SOURCE_DIR}/../include/annotation
                                            ${CMAKE_CURRENT_SOURCE_DIR}/../include/pipelines/detect
                                            ${CMAKE_CURRENT_SOURCE_DIR}/../include/pipelines/dataTypes
                                            ${VIENNA_RNA_INCLUDE_DIR}
                                            ${VIENNA_RNA_INCLUDE_DIR}/ViennaRNA
                                            ${HTSlib_INCLUDE_DIRS}
                                            ${HTSlib_INCLUDE_DIRS}/htslib)

target_link_libraries(tests     PRIVATE
                                GTest::gtest_main
                                seqan3::seqan3
                                Boost::program_options
                                Boost::math
                                Boost::uuid
                                segemehl)


include(GoogleTest)
gtest_discover_tests(tests)
