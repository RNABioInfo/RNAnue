# * Try to find ViennaRNA Once done, this will define
#
# VIENNA_RNA_FOUND - system has htslib VIENNA_RNA_INCLUDE_DIR - the htslib include directories
# VIENNA_RNA_LIBRARY - link these to use htslib

message(STATUS "Finding ViennaRNA")

find_path(
  VIENNA_RNA_INCLUDE_DIR
  NAMES ViennaRNA/cofold.h
  PATH_SUFFIXES include
)

find_library(
  VIENNA_RNA_LIBRARY
  NAMES libRNA.a
  PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu
)

if(VIENNA_RNA_INCLUDE_DIR AND VIENNA_RNA_LIBRARY)
    set(VIENNA_RNA_FOUND TRUE)
else()
    set(VIENNA_RNA_FOUND FALSE)
endif()

message(STATUS "   ViennaRNA include dirs: ${VIENNA_RNA_INCLUDE_DIR}")
message(STATUS "   ViennaRNA libraries: ${VIENNA_RNA_LIBRARY}")
message(STATUS "   ViennaRNA FOUND: ${VIENNA_RNA_FOUND} ")
