include(ExternalProject)

# ------------------------------------------------------------------------------
# Attempt to find system-installed HTSlib
# ------------------------------------------------------------------------------

find_package(HTSlib)
if(HTSlib_FOUND)
    message(STATUS "HTSlib_USE_STATIC_LIBS: ${HTSlib_USE_STATIC_LIBS}")
else()
    message(STATUS "HTSlib not found, building from source")

    # ------------------------------------------------------------------------------
    # Set up variables for building HTSlib
    # ------------------------------------------------------------------------------

    set(htslib_PREFIX ${CMAKE_BINARY_DIR}/submodules/htslib-prefix)
    set(htslib_INSTALL ${CMAKE_BINARY_DIR}/submodules/htslib-install)

    if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
        set(MAKE_COMMAND "$(MAKE)")
    else()
        find_program(MAKE_COMMAND NAMES make gmake)
    endif()

    message(STATUS "Building static htslib from source")
    message(NOTICE "Set ENV CFLAGS and CXXFLAGS if you use conda environment!")

    set(disable_flags --disable-gcs --disable-s3 --disable-plugins)

    # find lzma
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

    # config cmake files path deflate
    find_package(Deflate)
    if(Deflate_FOUND)
        include_directories(SYSTEM ${Deflate_INCLUDE_DIRS})
        list(APPEND deps_LIB ${Deflate_LIBRARIES})
    endif()

    message(STATUS " dependencies: ${deps_LIB}")
    message(STATUS " htslib make command: ${MAKE_COMMAND} $")

    message(STATUS "Configure command: ${CONFIGURE_COMMAND}")

    message(STATUS "ZLIB_BUILD: ${ZLIB_BUILD}")
    if(ZLIB_BUILD)
        include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/zlib.cmake)
        set(LOCAL_ZLIB_CONFIG "CPPFLAGS=-I${CMAKE_BINARY_DIR}/submodules/zlib-install/include/ LDFLAGS=-L${CMAKE_BINARY_DIR}/submodules/zlib-install/lib/")
        message(STATUS "Updated configure command: ${LOCAL_ZLIB_CONFIG}")
        list(APPEND deps_LIB ${zlib_LIBRARIES} ${CMAKE_BINARY_DIR}/submodules/zlib-install/lib/libz.a)
    else()
        find_package(ZLIB)
        if(ZLIB_FOUND)
            include_directories(SYSTEM ${ZLIB_INCLUDE_DIRS})
            list(APPEND deps_LIB ${ZLIB_LIBRARIES})
        else()
            set(ZLIB_BUILD TRUE)
            # build zlib from source
            message(STATUS "Building zlib from source")
            include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/zlib.cmake)
            set(LOCAL_ZLIB_CONFIG "CPPFLAGS=-I${CMAKE_BINARY_DIR}/submodules/zlib-install/include/ LDFLAGS=-L${CMAKE_BINARY_DIR}/submodules/zlib-install/lib/")
            message(STATUS "Updated configure command: ${LOCAL_ZLIB_CONFIG}")
            list(APPEND deps_LIB ${zlib_LIBRARIES} ${CMAKE_BINARY_DIR}/submodules/zlib-install/lib/libz.a)
        endif()
    endif()

    list(APPEND deps_LIB ${zlib_LIBRARIES})

    ExternalProject_Add(
        htslib
        PREFIX ${htslib_PREFIX}
        DOWNLOAD_EXTRACT_TIMESTAMP true
        URL https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2
        BUILD_IN_SOURCE 1
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND autoreconf -i && ./configure --prefix=${htslib_PREFIX} ${disable_flags} CXX=$ENV{CXX} CC=$ENV{CC} ${LOCAL_ZLIB_CONFIG}
        BUILD_COMMAND ${MAKE_COMMAND} lib-static CXX=$ENV{CXX} CC=$ENV{CC}
        INSTALL_COMMAND ${MAKE_COMMAND} install prefix=${htslib_INSTALL}
  )

    if (ZLIB_BUILD)
        add_dependencies(htslib zlib)
    endif()

    set(HTSlib_INCLUDE_DIRS ${htslib_INSTALL}/include ${htslib_INSTALL}/include/htslib ${CMAKE_BINARY_DIR}/submodules/zlib-install/include/)
    set(HTSlib_LIBRARIES ${htslib_INSTALL}/lib/libhts.a ${deps_LIB})

    message(STATUS "HTSlib_INCLUDE_DIRS: ${HTSlib_INCLUDE_DIRS}")
    message(STATUS "HTSlib_LIBRARIES: ${HTSlib_LIBRARIES}")

endif()
