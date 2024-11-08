cmake_minimum_required(VERSION 3.24)

project(RNAnue VERSION 2.0.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -fopenmp -Wno-psabi -Wno-interference-size")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Configure header file to pass the version number to the source code
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/Config.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/include/Config.hpp)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
# # # CMake Modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# # # SeqAn3
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan3/build_system)
find_package(seqan3 3.3.0 REQUIRED)

# # # Boost
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
set(Boost_DEBUG OFF)
set(Boost_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/submodules/boost)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/submodules/boost EXCLUDE_FROM_ALL)

# # # ViennaRNA
find_package(ViennaRNA)

# # # htslib
find_package(HTSlib)

# # # bam_sort
set(BAM_SORT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/bam_sort)

# # # segemehl
add_subdirectory(submodules/segemehl)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/parameters
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/annotation
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/pipelines/
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/pipelines/preprocess
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/pipelines/align
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/pipelines/detect
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/pipelines/analyze
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/pipelines/dataTypes
                    ${CMAKE_CURRENT_SOURCE_DIR}/submodules/bam_sort)

add_executable (RNAnue  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/Runner.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/annotation/FeatureAnnotator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/annotation/FeatureParser.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/annotation/FeatureWriter.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/annotation/Orientation.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/parameters/ParameterOptions.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/parameters/ParameterParser.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/dataTypes/GenomicRegion.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/preprocess/Preprocess.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/preprocess/PreprocessData.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/preprocess/Deduplicator.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/align/Align.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/align/AlignData.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/Detect.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/DetectData.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/SplitRecordsComplementarityEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/SplitRecordsEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/SplitRecordsHybridizationEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/SplitRecordsSplicingEvaluator.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/Analyze.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/AnalyzeData.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/RecordFragment.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/InteractionCluster.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/InteractionClusterGenerator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/ParallelInteractionClusterGenerator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/InteractionsWriter.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/StatisticEvaluator.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/Closing.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/CooptimalPairwiseAligner.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/Utility.cpp
)

target_sources(RNAnue PRIVATE
    ${BAM_SORT_DIR}/bam_sort.c
    ${BAM_SORT_DIR}/sam_opts.c
    ${BAM_SORT_DIR}/bedidx.c
    ${BAM_SORT_DIR}/bam.c
    ${BAM_SORT_DIR}/sam_utils.c
)

# # # htslib
if(NOT HTSlib_FOUND)
    message(STATUS "Building static htslib from source")
    include(htslib)
    add_dependencies(RNAnue htslib)
    add_dependencies(segemehl htslib)

    target_include_directories(segemehl PUBLIC include ${HTSlib_INCLUDE_DIRS})
else()
    message(STATUS "Using htslib from system")
endif()

target_link_libraries(RNAnue PRIVATE ${HTSlib_LIBRARIES})
target_include_directories(RNAnue PRIVATE ${HTSlib_INCLUDE_DIRS} ${HTSlib_INCLUDE_DIRS}/htslib)

# # # ViennaRNA
message(STATUS "Finding ViennaRNA status: ${VIENNA_RNA_FOUND}")
if (NOT VIENNA_RNA_FOUND)
    message(STATUS "Building static ViennaRNA from source")
    include(viennaRNA)
    add_dependencies(RNAnue viennaRNA)
else()
    message(STATUS "Using ViennaRNA from system")
endif()


target_link_libraries (RNAnue PRIVATE ${VIENNA_RNA_LIBRARY})
target_include_directories (RNAnue SYSTEM PRIVATE ${VIENNA_RNA_INCLUDE_DIR} ${VIENNA_RNA_INCLUDE_DIR}/ViennaRNA)

# # # Link Libraries
target_link_libraries (RNAnue PRIVATE   seqan3::seqan3
                                        Boost::program_options
                                        Boost::math
                                        Boost::uuid
                                        segemehl)

# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------
if (CMAKE_BUILD_TYPE STREQUAL DEBUG OR CMAKE_BUILD_TYPE STREQUAL Debug)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
    set(BUILD_TESTS ON)
    message(STATUS "Enabling tests")
    add_subdirectory(tests)

    if (NOT VIENNA_RNA_FOUND)
        message(STATUS "Adding ViennaRNA as dependency for tests")
        add_dependencies(tests viennaRNA)
    endif()

    if(NOT HTSlib_FOUND)
        message(STATUS "Adding htslib as dependency for tests")
        add_dependencies(tests htslib)
    endif()
endif()

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------
install(TARGETS RNAnue
    RUNTIME_DEPENDENCY_SET appDeps
)
install(RUNTIME_DEPENDENCY_SET appDeps
    PRE_EXCLUDE_REGEXES
        [[libc\.so\..*]]
        [[libgcc_s\.so\..*]]
        [[libm\.so\..*]]
        [[libstdc\+\+\.so\..*]]
)

# ----------------------------------------------------------------------------
# Packaging
# ----------------------------------------------------------------------------
set(CPACK_PACKAGE_NAME RNAnue)
set(CPACK_PACKAGE_VENDOR "RNA Biology & Bioinformatics; Institute of Biomedical Genetics; University of Stuttgart")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "RNAnue: efficient data analysis for RNA-RNA interactomics")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_VERBATIM_VARIABLES TRUE)
set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CPACK_PACKAGE_NAME})

set(CPACK_PACKAGE_CONTACT "Christopher Adelmann; RNA Biology & Bioinformatics; Institute of Biomedical Genetics; University of Stuttgart")
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE)

if(APPLE)
    message(STATUS "Building for MacOS")
    set(CPACK_GENERATOR TGZ productbuild)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(STATUS "Building for Linux")
    set(CPACK_GENERATOR TGZ DEB)
else()
    message(STATUS "Building for Other OS")
    set(CPACK_GENERATOR TGZ)
endif()

include(CPack)

cpack_add_component(RNAnue
    DISPLAY_NAME "RNAnue"
    DESCRIPTION "RNAnue: efficient data analysis for RNA-RNA interactomics"
    REQUIRED
)
