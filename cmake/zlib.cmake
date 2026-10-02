include(ExternalProject)

set(zlib_PREFIX ${CMAKE_BINARY_DIR}/submodules/zlib-prefix)
set(zlib_INSTALL ${CMAKE_BINARY_DIR}/submodules/zlib-install)
set(zlib_LIBRARIES ${zlib_INSTALL}/lib/libz.a)

ExternalProject_Add(
    zlib
    PREFIX ${zlib_PREFIX}
    DOWNLOAD_EXTRACT_TIMESTAMP true
    URL "https://github.com/madler/zlib/archive/refs/tags/v1.3.1.tar.gz"
    UPDATE_COMMAND ""
    # BUILD_IN_SOURCE 1 CONFIGURE_COMMAND ${zlib_PREFIX}/src/zlib/configure --prefix=${zlib_INSTALL} --static
    INSTALL_DIR ${zlib_INSTALL}
    CMAKE_ARGS
        "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
        "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${zlib_INSTALL} -DCMAKE_MACOSX_RPATH=0
    INSTALL_COMMAND ""
)

# Declare the installed archive at the step that actually creates it.
ExternalProject_Add_Step(zlib install_library
    COMMAND "${CMAKE_COMMAND}" --install <BINARY_DIR> --config "$<CONFIG>"
    DEPENDEES build
    DEPENDERS install
    BYPRODUCTS ${zlib_LIBRARIES}
)

include_directories(${zlib_INSTALL}/include)
