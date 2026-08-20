cmake_minimum_required(VERSION 3.24)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/select_compilers.cmake)
rnanue_select_compatible_compilers()

project(RNAnue VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/compiler_features.cmake)

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
include(dependency_providers)
rnanue_configure_dependency_provider(RNANUE_ZLIB_PROVIDER
    "zlib provider: AUTO tries installed zlib first, SYSTEM requires it, BUNDLED builds zlib from source when a bundled dependency needs it")
rnanue_detect_active_dependency_prefix(RNANUE_ACTIVE_DEPENDENCY_PREFIX)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux"
    AND RNANUE_ACTIVE_DEPENDENCY_PREFIX
    AND NOT RNANUE_DEPENDENCY_PROVIDER STREQUAL "BUNDLED")
  rnanue_prefer_active_prefix_libraries("${RNANUE_ACTIVE_DEPENDENCY_PREFIX}")
endif()

set(CPM_USE_LOCAL_PACKAGES ON CACHE BOOL "Try find_package before downloading CPM dependencies" FORCE)
if(RNANUE_DEPENDENCY_PROVIDER STREQUAL "SYSTEM")
  set(FETCHCONTENT_FULLY_DISCONNECTED ON CACHE BOOL
        "Disable FetchContent downloads in strict system dependency mode" FORCE)
  set(SEQAN3_NO_CEREAL ON CACHE BOOL
        "Disable optional SeqAn3 cereal support in strict system dependency mode" FORCE)
  find_path(SEQAN3_SDSL_INCLUDE_DIR
        NAMES sdsl/version.hpp
        PATH_SUFFIXES include
    )
  if(NOT SEQAN3_SDSL_INCLUDE_DIR)
    message(FATAL_ERROR
            "RNANUE_DEPENDENCY_PROVIDER=SYSTEM requires installed SDSL headers for SeqAn3. "
            "Install SDSL, set CMAKE_PREFIX_PATH/CMAKE_INCLUDE_PATH to its prefix, or use "
            "-DRNANUE_DEPENDENCY_PROVIDER=AUTO.")
  endif()
else()
  set(FETCHCONTENT_FULLY_DISCONNECTED OFF CACHE BOOL
        "Allow FetchContent downloads outside strict system dependency mode" FORCE)
endif()

# # # SeqAn3
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan3/cmake)
find_package(seqan3 3.3.0 REQUIRED)

# # # ViennaRNA
rnanue_configure_dependency_provider(RNANUE_VIENNARNA_PROVIDER
    "ViennaRNA provider: AUTO tries installed ViennaRNA first, SYSTEM requires it, BUNDLED builds ViennaRNA from source")
rnanue_provider_allows_system(RNANUE_VIENNARNA_PROVIDER RNANUE_VIENNARNA_ALLOW_SYSTEM)
rnanue_provider_allows_bundled(RNANUE_VIENNARNA_PROVIDER RNANUE_VIENNARNA_ALLOW_BUNDLED)
if(RNANUE_VIENNARNA_ALLOW_SYSTEM)
  find_package(ViennaRNA)
else()
  set(VIENNA_RNA_FOUND FALSE)
endif()
if(VIENNA_RNA_FOUND)
  rnanue_record_dependency("ViennaRNA" "${RNANUE_VIENNARNA_PROVIDER}" "${VIENNA_RNA_INCLUDE_DIR}" "${VIENNA_RNA_LIBRARIES}")
elseif(RNANUE_VIENNARNA_PROVIDER STREQUAL "SYSTEM")
  message(FATAL_ERROR
        "RNANUE_VIENNARNA_PROVIDER=SYSTEM requires an installed ViennaRNA package. "
        "Install ViennaRNA, set CMAKE_PREFIX_PATH to its prefix, or configure with "
        "-DRNANUE_VIENNARNA_PROVIDER=AUTO or BUNDLED.")
endif()

# # # htslib
rnanue_configure_dependency_provider(RNANUE_HTSLIB_PROVIDER
    "HTSlib provider: AUTO tries installed htslib first, SYSTEM requires it, BUNDLED builds htslib from source")
rnanue_provider_allows_system(RNANUE_HTSLIB_PROVIDER RNANUE_HTSLIB_ALLOW_SYSTEM)
rnanue_provider_allows_bundled(RNANUE_HTSLIB_PROVIDER RNANUE_HTSLIB_ALLOW_BUNDLED)
if(RNANUE_HTSLIB_ALLOW_SYSTEM)
  find_package(HTSlib)
else()
  set(HTSlib_FOUND FALSE)
endif()
if(HTSlib_FOUND)
  rnanue_record_dependency("HTSlib" "${RNANUE_HTSLIB_PROVIDER}" "${HTSlib_INCLUDE_DIRS}" "${HTSlib_LIBRARIES}")
elseif(RNANUE_HTSLIB_PROVIDER STREQUAL "SYSTEM")
  message(FATAL_ERROR
        "RNANUE_HTSLIB_PROVIDER=SYSTEM requires an installed htslib package. "
        "Install htslib, set CMAKE_PREFIX_PATH/HTSLIB_ROOT to its prefix, or configure with "
        "-DRNANUE_HTSLIB_PROVIDER=AUTO or BUNDLED.")
elseif(RNANUE_HTSLIB_ALLOW_BUNDLED)
  message(STATUS "Building static htslib from source")
  include(htslib)
  set(RNANUE_HTSLIB_EXTERNAL_TARGET htslib)
  rnanue_record_dependency("HTSlib" "${RNANUE_HTSLIB_PROVIDER}" "${HTSlib_INCLUDE_DIRS}" "${HTSlib_LIBRARIES}")
else()
  message(FATAL_ERROR "HTSlib was not found and bundled HTSlib is disabled.")
endif()

# # # Intel oneTBB
rnanue_configure_dependency_provider(RNANUE_TBB_PROVIDER
    "oneTBB provider: AUTO tries installed oneTBB first, SYSTEM requires it, BUNDLED builds oneTBB from source")
rnanue_provider_allows_system(RNANUE_TBB_PROVIDER RNANUE_TBB_ALLOW_SYSTEM)
rnanue_provider_allows_bundled(RNANUE_TBB_PROVIDER RNANUE_TBB_ALLOW_BUNDLED)
if(RNANUE_TBB_ALLOW_SYSTEM)
  rnanue_try_prefix_tbb(RNANUE_USED_PREFIX_TBB)
  if(NOT RNANUE_USED_PREFIX_TBB)
    find_package(TBB QUIET COMPONENTS tbb)
  endif()
else()
  set(TBB_FOUND FALSE)
endif()
if(TBB_FOUND)
  message(STATUS "Using oneTBB from installed package")
  rnanue_record_dependency("oneTBB" "${RNANUE_TBB_PROVIDER}" "" "TBB::tbb")
elseif(RNANUE_TBB_PROVIDER STREQUAL "SYSTEM")
  message(FATAL_ERROR
        "RNANUE_TBB_PROVIDER=SYSTEM requires installed oneTBB. Install oneTBB, set TBB_DIR, "
        "or configure with -DRNANUE_TBB_PROVIDER=AUTO or BUNDLED.")
elseif(RNANUE_TBB_ALLOW_BUNDLED)
  message(STATUS "Building bundled oneTBB from source")
  include(tbb)
endif()

# # # bam_sort
set(BAM_SORT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/submodules/bam_sort)

# # # segemehl
set(RNANUE_SAVED_DISABLE_FIND_PACKAGE_HTSlib "${CMAKE_DISABLE_FIND_PACKAGE_HTSlib}")
if(RNANUE_HTSLIB_EXTERNAL_TARGET)
  set(CMAKE_DISABLE_FIND_PACKAGE_HTSlib TRUE)
endif()
add_subdirectory(submodules/segemehl)
set(CMAKE_DISABLE_FIND_PACKAGE_HTSlib "${RNANUE_SAVED_DISABLE_FIND_PACKAGE_HTSlib}")

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
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/annotation/AnnotationHierarchyResolver.cpp
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
  add_dependencies(RNAnue ${RNANUE_HTSLIB_EXTERNAL_TARGET})
  add_dependencies(segemehl ${RNANUE_HTSLIB_EXTERNAL_TARGET})

  target_include_directories(segemehl PUBLIC include ${HTSlib_INCLUDE_DIRS})
else()
  message(STATUS "Using htslib from installed package")
endif()

target_link_libraries(RNAnue PRIVATE ${HTSlib_LIBRARIES})
target_include_directories(RNAnue PRIVATE ${HTSlib_INCLUDE_DIRS} ${HTSlib_INCLUDE_DIRS}/htslib)

# # # ViennaRNA
message(STATUS "Finding ViennaRNA status: ${VIENNA_RNA_FOUND}")
if (NOT VIENNA_RNA_FOUND)
  if(NOT RNANUE_VIENNARNA_ALLOW_BUNDLED)
    message(FATAL_ERROR "ViennaRNA was not found and bundled ViennaRNA is disabled.")
  endif()
  message(STATUS "Building static ViennaRNA from source")
  include(viennaRNA)
  add_dependencies(RNAnue viennaRNA)
  rnanue_record_dependency("ViennaRNA" "${RNANUE_VIENNARNA_PROVIDER}" "${VIENNA_RNA_INCLUDE_DIR}" "${VIENNA_RNA_LIBRARIES}")
else()
  message(STATUS "Using ViennaRNA from installed package")
endif()


# # # Boost
include(boost)
if(RNANUE_BOOST_EXTERNAL_TARGET)
  add_dependencies(RNAnue ${RNANUE_BOOST_EXTERNAL_TARGET})
endif()

# # # Matplotlib
include(matplotplusplus)

# # # Dependency Include Directories
target_include_directories (RNAnue SYSTEM PRIVATE   ${VIENNA_RNA_INCLUDE_DIR}
                                                    ${VIENNA_RNA_INCLUDE_DIR}/ViennaRNA
                                                    ${RNANUE_BOOST_INCLUDE_DIRS})

# # # Link Libraries
target_link_libraries (RNAnue PRIVATE   seqan3::seqan3
                                        segemehl
                                        TBB::tbb
                                        ${VIENNA_RNA_LIBRARIES}
                                        ${RNANUE_BOOST_LIBRARIES})

target_link_libraries (RNAnue PUBLIC matplot)

if(TARGET RNAnue_tbb)
  add_dependencies(RNAnue RNAnue_tbb)
endif()

rnanue_print_dependency_summary()

# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------
option(RNANUE_BUILD_TESTS "Build RNAnue tests" ON)

if (RNANUE_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)

  if (NOT VIENNA_RNA_FOUND)
    add_dependencies(tests viennaRNA)
  endif()

  if(NOT HTSlib_FOUND)
    add_dependencies(tests htslib)
  endif()

  if(RNANUE_BOOST_EXTERNAL_TARGET)
    add_dependencies(tests ${RNANUE_BOOST_EXTERNAL_TARGET})
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
