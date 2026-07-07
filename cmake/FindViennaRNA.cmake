# * Try to find ViennaRNA Once done, this will define
#
# VIENNA_RNA_FOUND - system has ViennaRNA
# VIENNA_RNA_INCLUDE_DIR - the ViennaRNA include directory
# VIENNA_RNA_LIBRARY - the ViennaRNA archive/library
# VIENNA_RNA_LIBRARIES - link these to use ViennaRNA

message(STATUS "Finding ViennaRNA")

find_path(
  VIENNA_RNA_INCLUDE_DIR
  NAMES ViennaRNA/cofold.h
  PATH_SUFFIXES include
)

find_library(
  VIENNA_RNA_LIBRARY
  NAMES RNA libRNA.a
  PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu
)

if(VIENNA_RNA_INCLUDE_DIR AND VIENNA_RNA_LIBRARY)
    set(VIENNA_RNA_FOUND TRUE)
    set(VIENNA_RNA_LIBRARIES ${VIENNA_RNA_LIBRARY})

    if(VIENNA_RNA_LIBRARY MATCHES "\\.(a|lib)$")
        set(VIENNA_RNA_REQUIRES_OPENMP FALSE)
        find_program(RNANUE_NM_TOOL NAMES ${CMAKE_NM} nm llvm-nm)

        if(RNANUE_NM_TOOL)
            execute_process(
                COMMAND ${RNANUE_NM_TOOL} -g "${VIENNA_RNA_LIBRARY}"
                OUTPUT_VARIABLE RNANUE_VIENNARNA_SYMBOLS
                ERROR_QUIET
            )

            if(RNANUE_VIENNARNA_SYMBOLS MATCHES "omp_set_dynamic|omp_get_|omp_set_|GOMP_")
                set(VIENNA_RNA_REQUIRES_OPENMP TRUE)
            endif()
        endif()

        if(VIENNA_RNA_REQUIRES_OPENMP)
            find_package(OpenMP QUIET COMPONENTS CXX)

            if(OpenMP_CXX_FOUND)
                list(APPEND VIENNA_RNA_LIBRARIES OpenMP::OpenMP_CXX)
                message(STATUS "   ViennaRNA static library requires OpenMP: linking OpenMP runtime")
            else()
                message(FATAL_ERROR
                    "ViennaRNA was found as static library '${VIENNA_RNA_LIBRARY}' and references "
                    "OpenMP symbols, but CMake could not find a C++ OpenMP runtime. Install the "
                    "matching OpenMP runtime for this compiler or configure with "
                    "-DRNANUE_VIENNARNA_PROVIDER=BUNDLED to build ViennaRNA without OpenMP.")
            endif()
        endif()
    endif()
else()
    set(VIENNA_RNA_FOUND FALSE)
    set(VIENNA_RNA_LIBRARIES "")
endif()

message(STATUS "   ViennaRNA include dirs: ${VIENNA_RNA_INCLUDE_DIR}")
message(STATUS "   ViennaRNA libraries: ${VIENNA_RNA_LIBRARY}")
message(STATUS "   ViennaRNA link libraries: ${VIENNA_RNA_LIBRARIES}")
message(STATUS "   ViennaRNA FOUND: ${VIENNA_RNA_FOUND} ")
