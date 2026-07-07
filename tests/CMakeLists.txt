# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------
include(FetchContent)
include(${CMAKE_CURRENT_LIST_DIR}/../cmake/dependency_providers.cmake)

# The top-level project still uses directory-wide include paths, but fetched
# dependencies should not inherit RNAnue's headers or Conda-injected SDK
# overlays. Keep this subdirectory explicit so GoogleTest is built in a clean
# compiler environment.
set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "")

foreach(RNANUE_TEST_FLAG_VAR CMAKE_C_FLAGS CMAKE_CXX_FLAGS)
    string(REGEX REPLACE
           "(^|[ \t]+)-isystem[ \t]+[^ \t]+/(miniconda3|miniforge3|anaconda3|mambaforge)/include($|[ \t]+)"
           " "
           ${RNANUE_TEST_FLAG_VAR}
           " ${${RNANUE_TEST_FLAG_VAR}} ")
    string(STRIP "${${RNANUE_TEST_FLAG_VAR}}" ${RNANUE_TEST_FLAG_VAR})
endforeach()

set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

rnanue_provider_allows_system(RNANUE_DEPENDENCY_PROVIDER RNANUE_GTEST_ALLOW_SYSTEM)
rnanue_provider_allows_bundled(RNANUE_DEPENDENCY_PROVIDER RNANUE_GTEST_ALLOW_BUNDLED)

if(RNANUE_GTEST_ALLOW_SYSTEM)
    find_package(GTest QUIET)
endif()

if(NOT TARGET GTest::gtest_main)
    if(RNANUE_DEPENDENCY_PROVIDER STREQUAL "SYSTEM")
        message(FATAL_ERROR
            "RNANUE_DEPENDENCY_PROVIDER=SYSTEM requires installed GoogleTest when "
            "RNANUE_BUILD_TESTS=ON. Install GTest, set CMAKE_PREFIX_PATH, disable tests "
            "with -DRNANUE_BUILD_TESTS=OFF, or use AUTO/BUNDLED dependency mode.")
    elseif(RNANUE_GTEST_ALLOW_BUNDLED)
        FetchContent_Declare(
            googletest
            DOWNLOAD_EXTRACT_TIMESTAMP true
            URL https://github.com/google/googletest/archive/1d17ea141d2c11b8917d2c7d029f1c4e2b9769b2.zip
        )

        FetchContent_MakeAvailable(googletest)
        rnanue_record_dependency("GoogleTest" "${RNANUE_DEPENDENCY_PROVIDER}" "" "GTest::gtest_main")
    endif()
else()
    rnanue_record_dependency("GoogleTest" "${RNANUE_DEPENDENCY_PROVIDER}" "" "GTest::gtest_main")
endif()

if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Newer Apple SDK Mach headers use the C11 spelling in headers that
    # GoogleTest includes from C++ sources. Clang accepts it; GCC needs the
    # C++ spelling. Scope the workaround to GoogleTest only.
    if(TARGET gtest)
        target_compile_definitions(gtest PRIVATE _Static_assert=static_assert)
    endif()
    if(TARGET gtest_main)
        target_compile_definitions(gtest_main PRIVATE _Static_assert=static_assert)
    endif()
endif()

add_executable(tests
    ${CMAKE_CURRENT_SOURCE_DIR}/ArmCoverageTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FeatureParserTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InteractionClusterTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InteractionClusterGeneratorTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InteractionsWriterTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ParallelInteractionClusterGeneratorTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SamFileUtilityTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/StatisticEvaluatorTests.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/FeatureAnnotator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/FeatureGrouper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/FeatureParser.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/annotation/GenomicFeatureTreeMerger.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/dataTypes/GenomicRegion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/dataTypes/GenomicOrientation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/dataTypes/RecordFragment.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/dataTypes/FeatureGroup.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/InteractionCluster.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/InteractionClusterGenerator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/InteractionsWriter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/ParallelInteractionClusterGenerator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/pipelines/analyze/StatisticEvaluator.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/../src/Utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/utility/SamFileUtility.cpp

    ${BAM_SORT_DIR}/bam_sort.c
    ${BAM_SORT_DIR}/sam_opts.c
    ${BAM_SORT_DIR}/bedidx.c
    ${BAM_SORT_DIR}/bam.c
    ${BAM_SORT_DIR}/sam_utils.c
)

target_include_directories(tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/annotation
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/concepts
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/parameters
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/parameters/options
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/parameters/options/optionTypes
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/pipelines/analyze
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/pipelines/dataTypes
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/utility
    ${BAM_SORT_DIR}
    ${VIENNA_RNA_INCLUDE_DIR}
    ${VIENNA_RNA_INCLUDE_DIR}/ViennaRNA
    ${HTSlib_INCLUDE_DIRS}
    ${HTSlib_INCLUDE_DIRS}/htslib
    ${RNANUE_BOOST_INCLUDE_DIRS}
)

target_link_libraries(tests PRIVATE
    GTest::gtest_main
    seqan3::seqan3
    segemehl
    TBB::tbb
    ${HTSlib_LIBRARIES}
    ${VIENNA_RNA_LIBRARIES}
    ${RNANUE_BOOST_LIBRARIES}
)

if(TARGET RNAnue_tbb)
    add_dependencies(tests RNAnue_tbb)
endif()

include(GoogleTest)
gtest_discover_tests(tests)
