include(CheckCXXSourceCompiles)
include(${CMAKE_CURRENT_LIST_DIR}/required_cxx_features.cmake)

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR
        "RNAnue requires GCC/G++ 14 or newer with C++23 standard-library support. "
        "The selected C++ compiler is '${CMAKE_CXX_COMPILER_ID}' (${CMAKE_CXX_COMPILER}). "
        "Reconfigure a fresh build directory with "
        "-DCMAKE_CXX_COMPILER=/path/to/g++-14 -DCMAKE_C_COMPILER=/path/to/gcc-14.")
endif()

if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14)
    message(FATAL_ERROR
        "RNAnue requires GCC/G++ 14 or newer. The selected compiler is "
        "GNU ${CMAKE_CXX_COMPILER_VERSION} (${CMAKE_CXX_COMPILER}). "
        "GCC 13 lacks required C++23 libstdc++ APIs used by RNAnue, including "
        "<print>, std::forward_like, std::ranges::to, and std::ranges::zip_view. "
        "Reconfigure a fresh build directory with "
        "-DCMAKE_CXX_COMPILER=/path/to/g++-14 -DCMAKE_C_COMPILER=/path/to/gcc-14.")
endif()

set(RNANUE_PREVIOUS_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++23")
check_cxx_source_compiles("${RNANUE_REQUIRED_CXX_FEATURES_SOURCE}"
                          RNANUE_HAS_REQUIRED_CXX_FEATURES)
set(CMAKE_REQUIRED_FLAGS "${RNANUE_PREVIOUS_REQUIRED_FLAGS}")
unset(RNANUE_PREVIOUS_REQUIRED_FLAGS)

if(NOT RNANUE_HAS_REQUIRED_CXX_FEATURES)
    message(FATAL_ERROR
        "The selected GCC/G++ compiler does not support the C++20/23 language and "
        "standard-library features required by RNAnue. Required features include "
        "std::format, std::formatter, std::print, std::forward_like, std::ranges::to, "
        "std::ranges::zip_view, std::unreachable, std::source_location, std::span, "
        "consteval, and coroutines. Use GCC/G++ 14 or newer and reconfigure a fresh "
        "build directory with -DCMAKE_CXX_COMPILER=/path/to/g++-14 "
        "-DCMAKE_C_COMPILER=/path/to/gcc-14.")
endif()

message(STATUS "RNAnue C++ compiler feature check passed")
