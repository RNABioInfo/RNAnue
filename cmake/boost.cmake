set(RNANUE_BOOST_PROVIDER "AUTO" CACHE STRING
    "Boost provider: AUTO tries system Boost first, SYSTEM requires it, BUNDLED builds Boost from source")
set_property(CACHE RNANUE_BOOST_PROVIDER PROPERTY STRINGS AUTO SYSTEM BUNDLED)

string(TOUPPER "${RNANUE_BOOST_PROVIDER}" RNANUE_BOOST_PROVIDER)
set(RNANUE_BOOST_PROVIDER "${RNANUE_BOOST_PROVIDER}" CACHE STRING
    "Boost provider: AUTO tries system Boost first, SYSTEM requires it, BUNDLED builds Boost from source" FORCE)

set(RNANUE_BOOST_PROVIDER_VALUES AUTO SYSTEM BUNDLED)
list(FIND RNANUE_BOOST_PROVIDER_VALUES "${RNANUE_BOOST_PROVIDER}" RNANUE_BOOST_PROVIDER_INDEX)
if(RNANUE_BOOST_PROVIDER_INDEX EQUAL -1)
    message(FATAL_ERROR
        "Invalid RNANUE_BOOST_PROVIDER='${RNANUE_BOOST_PROVIDER}'. "
        "Use AUTO, SYSTEM, or BUNDLED.")
endif()

set(RNANUE_BOOST_INCLUDE_DIRS "")
set(RNANUE_BOOST_LIBRARIES "")
set(RNANUE_BOOST_EXTERNAL_TARGET "")
set(RNANUE_CAN_USE_SYSTEM_BOOST TRUE)

if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(RNANUE_CAN_USE_SYSTEM_BOOST FALSE)
endif()

if(NOT RNANUE_BOOST_PROVIDER STREQUAL "BUNDLED" AND RNANUE_CAN_USE_SYSTEM_BOOST)
    find_package(Boost QUIET COMPONENTS program_options)

    if(Boost_FOUND)
        message(STATUS "Using Boost.Program_options from system")
        set(RNANUE_BOOST_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})

        if(TARGET Boost::program_options)
            set(RNANUE_BOOST_LIBRARIES Boost::program_options)
        else()
            set(RNANUE_BOOST_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
        endif()
    elseif(RNANUE_BOOST_PROVIDER STREQUAL "SYSTEM")
        message(FATAL_ERROR
            "RNANUE_BOOST_PROVIDER=SYSTEM requires Boost.Program_options. "
            "On Ubuntu install libboost-program-options-dev, or configure with "
            "-DRNANUE_BOOST_PROVIDER=AUTO or -DRNANUE_BOOST_PROVIDER=BUNDLED.")
    endif()
elseif(RNANUE_BOOST_PROVIDER STREQUAL "SYSTEM")
    message(FATAL_ERROR
        "RNANUE_BOOST_PROVIDER=SYSTEM is not supported with GCC on macOS because "
        "system Boost packages are commonly built with AppleClang/libc++ and are "
        "ABI-incompatible with RNAnue's GCC/libstdc++ build. Use "
        "-DRNANUE_BOOST_PROVIDER=BUNDLED on macOS, or use the release preset so "
        "AUTO can select the bundled fallback.")
elseif(NOT RNANUE_CAN_USE_SYSTEM_BOOST)
    message(STATUS "Skipping system Boost.Program_options because GCC on macOS requires a compatible Boost build")
endif()

if(NOT RNANUE_BOOST_LIBRARIES)
    include(ExternalProject)

    set(RNANUE_BUNDLED_BOOST_ROOT ${CMAKE_BINARY_DIR}/submodules/boost-prefix)
    set(RNANUE_BUNDLED_BOOST_INSTALL ${CMAKE_BINARY_DIR}/submodules/boost-install)
    set(RNANUE_BUNDLED_BOOST_INCLUDE_DIR ${RNANUE_BUNDLED_BOOST_INSTALL}/include)
    set(RNANUE_BUNDLED_BOOST_LIB_DIR ${RNANUE_BUNDLED_BOOST_INSTALL}/lib)
    set(RNANUE_BUNDLED_BOOST_TOOLSET gcc)
    set(RNANUE_BUNDLED_BOOST_CXXFLAGS "")

    if(APPLE)
        if(CMAKE_OSX_SYSROOT)
            list(APPEND RNANUE_BUNDLED_BOOST_CXXFLAGS "-isysroot" "${CMAKE_OSX_SYSROOT}")
        endif()

        if(CMAKE_OSX_DEPLOYMENT_TARGET)
            list(APPEND RNANUE_BUNDLED_BOOST_CXXFLAGS "-mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}")
        endif()

        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            list(APPEND RNANUE_BUNDLED_BOOST_CXXFLAGS "-D_Static_assert=static_assert")
        endif()
    endif()

    string(JOIN " " RNANUE_BUNDLED_BOOST_CXXFLAGS_STRING ${RNANUE_BUNDLED_BOOST_CXXFLAGS})

    message(STATUS "Building bundled Boost.Program_options from source")
    message(STATUS "Using Boost cxx compiler: ${CMAKE_CXX_COMPILER}")
    if(RNANUE_BUNDLED_BOOST_CXXFLAGS_STRING)
        message(STATUS "Using Boost cxx flags: ${RNANUE_BUNDLED_BOOST_CXXFLAGS_STRING}")
    endif()

    ExternalProject_Add(
        Boost
        PREFIX ${RNANUE_BUNDLED_BOOST_ROOT}
        BUILD_IN_SOURCE 1
        DOWNLOAD_EXTRACT_TIMESTAMP true
        URL "https://github.com/boostorg/boost/releases/download/boost-1.86.0/boost-1.86.0-cmake.tar.gz"
        CONFIGURE_COMMAND ${CMAKE_COMMAND}
            "-DRNANUE_BOOST_STEP=bootstrap"
            "-DRNANUE_BOOST_SOURCE_DIR=<SOURCE_DIR>"
            "-DRNANUE_BOOST_INSTALL_DIR=<INSTALL_DIR>"
            "-DRNANUE_BOOST_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
            "-DRNANUE_BOOST_CXXFLAGS=${RNANUE_BUNDLED_BOOST_CXXFLAGS_STRING}"
            "-DRNANUE_BOOST_TOOLSET=${RNANUE_BUNDLED_BOOST_TOOLSET}"
            -P "${CMAKE_CURRENT_LIST_DIR}/boost_build_step.cmake"
        BUILD_COMMAND ${CMAKE_COMMAND}
            "-DRNANUE_BOOST_STEP=build"
            "-DRNANUE_BOOST_SOURCE_DIR=<SOURCE_DIR>"
            "-DRNANUE_BOOST_INSTALL_DIR=<INSTALL_DIR>"
            "-DRNANUE_BOOST_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
            "-DRNANUE_BOOST_CXXFLAGS=${RNANUE_BUNDLED_BOOST_CXXFLAGS_STRING}"
            "-DRNANUE_BOOST_TOOLSET=${RNANUE_BUNDLED_BOOST_TOOLSET}"
            -P "${CMAKE_CURRENT_LIST_DIR}/boost_build_step.cmake"
        INSTALL_COMMAND ""
        INSTALL_DIR ${RNANUE_BUNDLED_BOOST_INSTALL}
    )

    set(RNANUE_BOOST_INCLUDE_DIRS ${RNANUE_BUNDLED_BOOST_INCLUDE_DIR})
    set(RNANUE_BOOST_LIBRARIES ${RNANUE_BUNDLED_BOOST_LIB_DIR}/libboost_program_options.a)
    set(RNANUE_BOOST_EXTERNAL_TARGET Boost)

    message(STATUS "Bundled Boost static libs: ${RNANUE_BOOST_LIBRARIES}")
    message(STATUS "Bundled Boost include dir: ${RNANUE_BOOST_INCLUDE_DIRS}")
endif()
