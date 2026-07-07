include(ExternalProject)
include(${CMAKE_CURRENT_LIST_DIR}/dependency_providers.cmake)

# ------------------------------------------------------------------------------
# Build bundled HTSlib
# ------------------------------------------------------------------------------

set(htslib_PREFIX ${CMAKE_BINARY_DIR}/submodules/htslib-prefix)
set(htslib_INSTALL ${CMAKE_BINARY_DIR}/submodules/htslib-install)
set(deps_LIB "")
set(LOCAL_ZLIB_CONFIG "")
set(ZLIB_BUILD FALSE)

if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(MAKE_COMMAND "$(MAKE)")
else()
    find_program(MAKE_COMMAND NAMES make gmake)
endif()

message(STATUS "Building static htslib from source")

set(disable_flags --disable-gcs --disable-s3 --disable-plugins)

find_package(LibLZMA)
if(LIBLZMA_FOUND)
    include_directories(SYSTEM ${LIBLZMA_INCLUDE_DIRS})
    list(APPEND deps_LIB ${LIBLZMA_LIBRARIES})
else()
    list(APPEND disable_flags --disable-lzma)
endif()

find_package(CURL)
if(CURL_FOUND)
    include_directories(SYSTEM ${CURL_INCLUDE_DIRS})
    list(APPEND deps_LIB ${CURL_LIBRARIES})
else()
    list(APPEND disable_flags --disable-libcurl)
endif()

find_package(BZip2)
if(BZIP2_FOUND)
    include_directories(SYSTEM ${BZIP2_INCLUDE_DIRS})
    list(APPEND deps_LIB ${BZIP2_LIBRARIES})
else()
    list(APPEND disable_flags --disable-bz2)
endif()

find_package(Deflate)
if(Deflate_FOUND)
    include_directories(SYSTEM ${Deflate_INCLUDE_DIRS})
    list(APPEND deps_LIB ${Deflate_LIBRARIES})
endif()

rnanue_configure_dependency_provider(RNANUE_ZLIB_PROVIDER
    "zlib provider: AUTO tries installed zlib first, SYSTEM requires it, BUNDLED builds zlib from source")
rnanue_provider_allows_system(RNANUE_ZLIB_PROVIDER RNANUE_ZLIB_ALLOW_SYSTEM)
rnanue_provider_allows_bundled(RNANUE_ZLIB_PROVIDER RNANUE_ZLIB_ALLOW_BUNDLED)

if(RNANUE_ZLIB_ALLOW_SYSTEM)
    find_package(ZLIB)
endif()

if(ZLIB_FOUND)
    include_directories(SYSTEM ${ZLIB_INCLUDE_DIRS})
    list(APPEND deps_LIB ${ZLIB_LIBRARIES})
    rnanue_record_dependency("zlib" "${RNANUE_ZLIB_PROVIDER}" "${ZLIB_INCLUDE_DIRS}" "${ZLIB_LIBRARIES}")
elseif(RNANUE_ZLIB_PROVIDER STREQUAL "SYSTEM")
    message(FATAL_ERROR
        "RNANUE_ZLIB_PROVIDER=SYSTEM requires installed zlib. Install zlib, set ZLIB_ROOT, "
        "or configure with -DRNANUE_ZLIB_PROVIDER=AUTO or BUNDLED.")
elseif(RNANUE_ZLIB_ALLOW_BUNDLED)
    set(ZLIB_BUILD TRUE)
    message(STATUS "Building bundled zlib from source")
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/zlib.cmake)
    set(LOCAL_ZLIB_CONFIG "CPPFLAGS=-I${zlib_INSTALL}/include/" "LDFLAGS=-L${zlib_INSTALL}/lib/")
    list(APPEND deps_LIB ${zlib_LIBRARIES})
    rnanue_record_dependency("zlib" "${RNANUE_ZLIB_PROVIDER}" "${zlib_INSTALL}/include" "${zlib_LIBRARIES}")
else()
    message(FATAL_ERROR "zlib was not found and bundled zlib is disabled.")
endif()

message(STATUS "HTSlib dependencies: ${deps_LIB}")
message(STATUS "HTSlib make command: ${MAKE_COMMAND}")

ExternalProject_Add(
    htslib
    PREFIX ${htslib_PREFIX}
    DOWNLOAD_EXTRACT_TIMESTAMP true
    URL https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND autoreconf -i && ./configure --prefix=${htslib_PREFIX} ${disable_flags} "CXX=${CMAKE_CXX_COMPILER}" "CC=${CMAKE_C_COMPILER}" ${LOCAL_ZLIB_CONFIG}
    BUILD_COMMAND ${MAKE_COMMAND} "CXX=${CMAKE_CXX_COMPILER}" "CC=${CMAKE_C_COMPILER}" lib-static
    INSTALL_COMMAND ${MAKE_COMMAND} install "CXX=${CMAKE_CXX_COMPILER}" "CC=${CMAKE_C_COMPILER}" prefix=${htslib_INSTALL}
)

if(ZLIB_BUILD)
    add_dependencies(htslib zlib)
endif()

set(HTSlib_INCLUDE_DIRS ${htslib_INSTALL}/include ${htslib_INSTALL}/include/htslib)
if(ZLIB_BUILD)
    list(APPEND HTSlib_INCLUDE_DIRS ${zlib_INSTALL}/include)
endif()
set(HTSlib_LIBRARIES ${htslib_INSTALL}/lib/libhts.a ${deps_LIB})

message(STATUS "HTSlib_INCLUDE_DIRS: ${HTSlib_INCLUDE_DIRS}")
message(STATUS "HTSlib_LIBRARIES: ${HTSlib_LIBRARIES}")
