include(ExternalProject)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_DEBUG OFF)

set(Boost_ROOT ${CMAKE_BINARY_DIR}/submodules/boost-prefix)

set(Boost_INSTALL ${CMAKE_BINARY_DIR}/submodules/boost-install )
set(Boost_INCLUDE_DIR ${Boost_INSTALL}/include )
set(Boost_LIB_DIR ${Boost_INSTALL}/lib )

ExternalProject_Add(
      Boost
      PREFIX ${Boost_ROOT}
      BUILD_IN_SOURCE 1
      DOWNLOAD_EXTRACT_TIMESTAMP true
      URL "https://github.com/boostorg/boost/releases/download/boost-1.86.0/boost-1.86.0-cmake.tar.gz"
      CONFIGURE_COMMAND ./bootstrap.sh
      --prefix=<INSTALL_DIR>
      --with-toolset=gcc
      --with-libraries=program_options
      --with-libraries=math
      BUILD_COMMAND ./b2 install --cxx=${CMAKE_CXX_COMPILER}  link=static variant=release threading=multi runtime-link=static
      INSTALL_COMMAND ""
      INSTALL_DIR ${Boost_INSTALL}
)

set(Boost_LIBRARIES
    ${Boost_LIB_DIR}/libboost_program_options.a
    ${Boost_LIB_DIR}/libboost_math_c99.a
    ${Boost_LIB_DIR}/libboost_math_c99f.a
    ${Boost_LIB_DIR}/libboost_math_c99l.a
    ${Boost_LIB_DIR}/libboost_math_tr1.a
    ${Boost_LIB_DIR}/libboost_math_tr1f.a
    ${Boost_LIB_DIR}/libboost_math_tr1l.a
)

message( STATUS "Boost static libs: " ${Boost_LIBRARIES} )
message( STATUS "Boost include dir: "  ${Boost_INCLUDE_DIR} )
