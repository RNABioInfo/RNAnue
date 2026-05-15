cmake_minimum_required(VERSION 3.24)

if(CMAKE_HOST_APPLE AND NOT DEFINED CMAKE_CXX_COMPILER)
    if(DEFINED ENV{CXX} AND NOT "$ENV{CXX}" MATCHES "^/usr/bin/"
       AND "$ENV{CXX}" MATCHES "g\\+\\+(-[0-9]+)?$")
        if(IS_ABSOLUTE "$ENV{CXX}")
            set(RNANUE_MACOS_GXX "$ENV{CXX}")
        else()
            find_program(RNANUE_MACOS_GXX
                         NAMES "$ENV{CXX}"
                         HINTS /opt/homebrew/bin /usr/local/bin /opt/local/bin)
        endif()
    else()
        find_program(RNANUE_MACOS_GXX
                     NAMES g++-14 g++-13 g++-12
                     HINTS /opt/homebrew/bin /usr/local/bin /opt/local/bin)
    endif()

    if(RNANUE_MACOS_GXX)
        set(CMAKE_CXX_COMPILER "${RNANUE_MACOS_GXX}" CACHE FILEPATH "C++ compiler" FORCE)

        get_filename_component(RNANUE_MACOS_GXX_DIR "${RNANUE_MACOS_GXX}" DIRECTORY)
        get_filename_component(RNANUE_MACOS_GXX_NAME "${RNANUE_MACOS_GXX}" NAME)
        string(REPLACE "g++" "gcc" RNANUE_MACOS_GCC_NAME "${RNANUE_MACOS_GXX_NAME}")
        find_program(RNANUE_MACOS_GCC
                     NAMES "${RNANUE_MACOS_GCC_NAME}"
                     HINTS "${RNANUE_MACOS_GXX_DIR}"
                     NO_DEFAULT_PATH)

        if(RNANUE_MACOS_GCC)
            set(CMAKE_C_COMPILER "${RNANUE_MACOS_GCC}" CACHE FILEPATH "C compiler" FORCE)
        endif()
    endif()
endif()

project(RNAnue VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR
        "RNAnue requires GCC with C++23 support on macOS. Reconfigure a fresh build directory "
        "with -DCMAKE_CXX_COMPILER=/path/to/g++ and -DCMAKE_C_COMPILER=/path/to/gcc, "
        "or set CXX and CC before the first CMake configure step.")
endif()

if(APPLE AND NOT CMAKE_C_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR
        "RNAnue requires a GCC C compiler on macOS. Reconfigure a fresh build directory "
        "with -DCMAKE_C_COMPILER=/path/to/gcc, or set CC before the first CMake configure step.")
endif()

if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    foreach(RNANUE_CXX_FLAG_VAR
            CMAKE_CXX_FLAGS
            CMAKE_CXX_FLAGS_DEBUG
            CMAKE_CXX_FLAGS_RELEASE
            CMAKE_CXX_FLAGS_RELWITHDEBINFO
            CMAKE_CXX_FLAGS_MINSIZEREL)
        string(REPLACE "-stdlib=libc++" "" ${RNANUE_CXX_FLAG_VAR} "${${RNANUE_CXX_FLAG_VAR}}")
        set(${RNANUE_CXX_FLAG_VAR} "${${RNANUE_CXX_FLAG_VAR}}" CACHE STRING "" FORCE)
    endforeach()
endif()

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(CMAKE_CXX_COMPILER MATCHES "^(/opt/homebrew|/usr/local|/opt/local)/")
    get_filename_component(RNANUE_COMPILER_BIN_DIR "${CMAKE_CXX_COMPILER}" DIRECTORY)
    get_filename_component(RNANUE_COMPILER_PREFIX "${RNANUE_COMPILER_BIN_DIR}" DIRECTORY)
    list(PREPEND CMAKE_PREFIX_PATH "${RNANUE_COMPILER_PREFIX}")
endif()

# Configure header file to pass the version number to the source code
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/Config.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/include/Config.hpp)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
# # # CMake Modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# # # SeqAn3
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan3/cmake)
find_package(seqan3 3.3.0 REQUIRED)

# # # ViennaRNA
find_package(ViennaRNA)

# # # htslib
find_package(HTSlib)

# # # Intel oneTBB
find_package(TBB REQUIRED COMPONENTS tbb)

# # # bam_sort
set(BAM_SORT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/bam_sort)

# # # segemehl
add_subdirectory(submodules/segemehl)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/concepts
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/utility
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/plotting
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/parameters
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/parameters/options
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/parameters/options/optionTypes
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/annotation
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/pipelines
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/pipelines/preprocess
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/pipelines/align
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/pipelines/detect
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/pipelines/analyze
                    ${CMAKE_CURRENT_SOURCE_DIR}/include/pipelines/dataTypes
                    ${BAM_SORT_DIR})

add_executable (RNAnue  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/Runner.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/annotation/FeatureAnnotator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/annotation/FeatureGrouper.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/annotation/FeatureParser.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/annotation/FeatureWriter.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/annotation/GenomicFeatureTreeMerger.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/utility/SamFileUtility.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/utility/SequenceFileUtility.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/parameters/ParameterOptions.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/parameters/ParameterParser.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/dataTypes/GenomicRegion.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/dataTypes/RecordFragment.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/dataTypes/GenomicOrientation.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/dataTypes/FeatureGroup.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/preprocess/Preprocess.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/preprocess/PreprocessData.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/preprocess/Deduplicator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/preprocess/DeduplicatorBySequenceSingleEnd.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/preprocess/DeduplicatorBySequencePairedEnd.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/preprocess/PairedEndPreprocessor.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/preprocess/SingleEndPreprocessor.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/align/Align.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/align/AlignData.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/align/GeneCopyMasker.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/CrosslinkingSitesEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/Detect.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/DetectData.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/ComplementarityEvaluationStep.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/HitGroupEvaluator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/HybridizationEvaluationStep.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/detect/SplicingEvaluationStep.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/Analyze.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/AnalyzeData.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/InteractionCluster.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/InteractionClusterGenerator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/ParallelInteractionClusterGenerator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/InteractionsWriter.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/pipelines/analyze/StatisticEvaluator.cpp

                        ${CMAKE_CURRENT_SOURCE_DIR}/src/Closing.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/Utility.cpp
)

target_compile_options(RNAnue PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall>
    $<$<COMPILE_LANGUAGE:CXX>:-Wextra>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-psabi>
    $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-interference-size>
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


# # # Boost
include(boost)
add_dependencies(RNAnue Boost)

# # # Matplotlib
include(matplotplusplus)
add_dependencies(RNAnue matplot)

# # # Dependency Include Directories
target_include_directories (RNAnue SYSTEM PRIVATE   ${VIENNA_RNA_INCLUDE_DIR}
                                                    ${VIENNA_RNA_INCLUDE_DIR}/ViennaRNA
                                                    ${Boost_INCLUDE_DIR})

# # # Link Libraries
target_link_libraries (RNAnue PRIVATE   seqan3::seqan3
                                        segemehl
                                        TBB::tbb
                                        ${VIENNA_RNA_LIBRARY}
                                        ${Boost_LIBRARIES})

target_link_libraries (RNAnue PUBLIC matplot)

# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------
if (CMAKE_BUILD_TYPE STREQUAL DEBUG OR CMAKE_BUILD_TYPE STREQUAL Debug)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
    # set(BUILD_TESTS ON)

    # message(STATUS "--------------")
    # message(STATUS "Enabling tests")
    # message(STATUS "--------------")

    # add_subdirectory(tests)

    # if (NOT VIENNA_RNA_FOUND)
    #     message(STATUS "Adding ViennaRNA as dependency for tests")
    #     add_dependencies(tests viennaRNA)
    # endif()

    # if(NOT HTSlib_FOUND)
    #     message(STATUS "Adding htslib as dependency for tests")
    #     add_dependencies(tests htslib)
    # endif()

    # add_dependencies(tests Boost)
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
