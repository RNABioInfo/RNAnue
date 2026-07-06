include_guard(GLOBAL)

include("${CMAKE_CURRENT_LIST_DIR}/required_cxx_features.cmake")

set(RNANUE_MINIMUM_GNU_COMPILER_VERSION 14)

function(rnanue_resolve_compiler RNANUE_COMPILER_CANDIDATE RNANUE_OUT)
    if(NOT RNANUE_COMPILER_CANDIDATE)
        set(${RNANUE_OUT} "" PARENT_SCOPE)
        return()
    endif()

    if(IS_ABSOLUTE "${RNANUE_COMPILER_CANDIDATE}")
        if(EXISTS "${RNANUE_COMPILER_CANDIDATE}")
            set(${RNANUE_OUT} "${RNANUE_COMPILER_CANDIDATE}" PARENT_SCOPE)
        else()
            set(${RNANUE_OUT} "" PARENT_SCOPE)
        endif()
        return()
    endif()

    find_program(RNANUE_RESOLVED_COMPILER
                 NAMES "${RNANUE_COMPILER_CANDIDATE}"
                 NO_CACHE)
    set(${RNANUE_OUT} "${RNANUE_RESOLVED_COMPILER}" PARENT_SCOPE)
endfunction()

function(rnanue_probe_cxx_compiler
         RNANUE_CXX_CANDIDATE
         RNANUE_OK_OUT
         RNANUE_RESOLVED_OUT
         RNANUE_VERSION_OUT
         RNANUE_REASON_OUT)
    rnanue_resolve_compiler("${RNANUE_CXX_CANDIDATE}" RNANUE_RESOLVED_CXX)

    if(NOT RNANUE_RESOLVED_CXX)
        set(${RNANUE_OK_OUT} FALSE PARENT_SCOPE)
        set(${RNANUE_RESOLVED_OUT} "" PARENT_SCOPE)
        set(${RNANUE_VERSION_OUT} "" PARENT_SCOPE)
        set(${RNANUE_REASON_OUT}
            "compiler '${RNANUE_CXX_CANDIDATE}' was not found"
            PARENT_SCOPE)
        return()
    endif()

    execute_process(
        COMMAND "${RNANUE_RESOLVED_CXX}" --version
        RESULT_VARIABLE RNANUE_VERSION_RESULT
        OUTPUT_VARIABLE RNANUE_VERSION_BANNER
        ERROR_VARIABLE RNANUE_VERSION_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE)

    if(NOT RNANUE_VERSION_RESULT EQUAL 0)
        set(${RNANUE_OK_OUT} FALSE PARENT_SCOPE)
        set(${RNANUE_RESOLVED_OUT} "${RNANUE_RESOLVED_CXX}" PARENT_SCOPE)
        set(${RNANUE_VERSION_OUT} "" PARENT_SCOPE)
        set(${RNANUE_REASON_OUT}
            "could not execute '${RNANUE_RESOLVED_CXX} --version': ${RNANUE_VERSION_ERROR}"
            PARENT_SCOPE)
        return()
    endif()

    if(RNANUE_VERSION_BANNER MATCHES "[Cc]lang")
        set(${RNANUE_OK_OUT} FALSE PARENT_SCOPE)
        set(${RNANUE_RESOLVED_OUT} "${RNANUE_RESOLVED_CXX}" PARENT_SCOPE)
        set(${RNANUE_VERSION_OUT} "" PARENT_SCOPE)
        set(${RNANUE_REASON_OUT}
            "compiler is Clang/AppleClang, but RNAnue requires GCC/G++ ${RNANUE_MINIMUM_GNU_COMPILER_VERSION} or newer"
            PARENT_SCOPE)
        return()
    endif()

    execute_process(
        COMMAND "${RNANUE_RESOLVED_CXX}" -dumpfullversion -dumpversion
        RESULT_VARIABLE RNANUE_DUMP_VERSION_RESULT
        OUTPUT_VARIABLE RNANUE_DUMP_VERSION
        ERROR_VARIABLE RNANUE_DUMP_VERSION_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE)

    string(REGEX MATCH "[0-9]+(\\.[0-9]+)*" RNANUE_GNU_VERSION "${RNANUE_DUMP_VERSION}")
    if(NOT RNANUE_GNU_VERSION)
        string(REGEX MATCH "[0-9]+(\\.[0-9]+)*" RNANUE_GNU_VERSION "${RNANUE_VERSION_BANNER}")
    endif()

    if(NOT RNANUE_DUMP_VERSION_RESULT EQUAL 0 OR NOT RNANUE_GNU_VERSION)
        set(${RNANUE_OK_OUT} FALSE PARENT_SCOPE)
        set(${RNANUE_RESOLVED_OUT} "${RNANUE_RESOLVED_CXX}" PARENT_SCOPE)
        set(${RNANUE_VERSION_OUT} "" PARENT_SCOPE)
        set(${RNANUE_REASON_OUT}
            "could not determine a GNU compiler version: ${RNANUE_DUMP_VERSION_ERROR}"
            PARENT_SCOPE)
        return()
    endif()

    if(RNANUE_GNU_VERSION VERSION_LESS "${RNANUE_MINIMUM_GNU_COMPILER_VERSION}")
        set(${RNANUE_OK_OUT} FALSE PARENT_SCOPE)
        set(${RNANUE_RESOLVED_OUT} "${RNANUE_RESOLVED_CXX}" PARENT_SCOPE)
        set(${RNANUE_VERSION_OUT} "${RNANUE_GNU_VERSION}" PARENT_SCOPE)
        set(${RNANUE_REASON_OUT}
            "GNU ${RNANUE_GNU_VERSION} is older than required GCC/G++ ${RNANUE_MINIMUM_GNU_COMPILER_VERSION}"
            PARENT_SCOPE)
        return()
    endif()

    set(RNANUE_PROBE_DIR "${CMAKE_BINARY_DIR}/CMakeFiles/rnanue-compiler-probe")
    set(RNANUE_PROBE_SOURCE "${RNANUE_PROBE_DIR}/required_cxx_features.cpp")
    set(RNANUE_PROBE_BINARY "${RNANUE_PROBE_DIR}/required_cxx_features")
    file(MAKE_DIRECTORY "${RNANUE_PROBE_DIR}")
    file(WRITE "${RNANUE_PROBE_SOURCE}" "${RNANUE_REQUIRED_CXX_FEATURES_SOURCE}")
    file(REMOVE "${RNANUE_PROBE_BINARY}")

    execute_process(
        COMMAND "${RNANUE_RESOLVED_CXX}" -std=c++23 "${RNANUE_PROBE_SOURCE}" -o "${RNANUE_PROBE_BINARY}"
        RESULT_VARIABLE RNANUE_PROBE_RESULT
        OUTPUT_VARIABLE RNANUE_PROBE_OUTPUT
        ERROR_VARIABLE RNANUE_PROBE_ERROR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE)

    if(NOT RNANUE_PROBE_RESULT EQUAL 0)
        if(RNANUE_PROBE_ERROR)
            string(REGEX REPLACE "[\r\n].*" "" RNANUE_PROBE_FIRST_ERROR "${RNANUE_PROBE_ERROR}")
        else()
            string(REGEX REPLACE "[\r\n].*" "" RNANUE_PROBE_FIRST_ERROR "${RNANUE_PROBE_OUTPUT}")
        endif()

        set(${RNANUE_OK_OUT} FALSE PARENT_SCOPE)
        set(${RNANUE_RESOLVED_OUT} "${RNANUE_RESOLVED_CXX}" PARENT_SCOPE)
        set(${RNANUE_VERSION_OUT} "${RNANUE_GNU_VERSION}" PARENT_SCOPE)
        set(${RNANUE_REASON_OUT}
            "failed RNAnue's C++23 feature probe: ${RNANUE_PROBE_FIRST_ERROR}"
            PARENT_SCOPE)
        return()
    endif()

    set(${RNANUE_OK_OUT} TRUE PARENT_SCOPE)
    set(${RNANUE_RESOLVED_OUT} "${RNANUE_RESOLVED_CXX}" PARENT_SCOPE)
    set(${RNANUE_VERSION_OUT} "${RNANUE_GNU_VERSION}" PARENT_SCOPE)
    set(${RNANUE_REASON_OUT} "compatible" PARENT_SCOPE)
endfunction()

function(rnanue_sanitize_libcxx_flag_from_string RNANUE_INPUT RNANUE_OUTPUT)
    string(REPLACE "-stdlib=libc++" "" RNANUE_SANITIZED "${RNANUE_INPUT}")
    string(REGEX REPLACE "  +" " " RNANUE_SANITIZED "${RNANUE_SANITIZED}")
    string(STRIP "${RNANUE_SANITIZED}" RNANUE_SANITIZED)
    set(${RNANUE_OUTPUT} "${RNANUE_SANITIZED}" PARENT_SCOPE)
endfunction()

function(rnanue_sanitize_gnu_compiler_flags)
    foreach(RNANUE_ENV_VAR CFLAGS CXXFLAGS)
        if(DEFINED ENV{${RNANUE_ENV_VAR}} AND NOT "$ENV{${RNANUE_ENV_VAR}}" STREQUAL "")
            set(RNANUE_ENV_FLAGS "$ENV{${RNANUE_ENV_VAR}}")
            if(RNANUE_ENV_FLAGS MATCHES "-stdlib=libc\\+\\+")
                rnanue_sanitize_libcxx_flag_from_string("${RNANUE_ENV_FLAGS}" RNANUE_SANITIZED_FLAGS)
                set(ENV{${RNANUE_ENV_VAR}} "${RNANUE_SANITIZED_FLAGS}")
                message(STATUS
                    "Removed -stdlib=libc++ from ${RNANUE_ENV_VAR} for GCC/libstdc++ compatibility")
            endif()
        endif()
    endforeach()

    foreach(RNANUE_FLAG_VAR
            CMAKE_C_FLAGS
            CMAKE_CXX_FLAGS
            CMAKE_C_FLAGS_DEBUG
            CMAKE_CXX_FLAGS_DEBUG
            CMAKE_C_FLAGS_RELEASE
            CMAKE_CXX_FLAGS_RELEASE
            CMAKE_C_FLAGS_RELWITHDEBINFO
            CMAKE_CXX_FLAGS_RELWITHDEBINFO
            CMAKE_C_FLAGS_MINSIZEREL
            CMAKE_CXX_FLAGS_MINSIZEREL)
        if(DEFINED ${RNANUE_FLAG_VAR} AND NOT "${${RNANUE_FLAG_VAR}}" STREQUAL ""
           AND "${${RNANUE_FLAG_VAR}}" MATCHES "-stdlib=libc\\+\\+")
            rnanue_sanitize_libcxx_flag_from_string("${${RNANUE_FLAG_VAR}}" RNANUE_SANITIZED_FLAGS)
            set(${RNANUE_FLAG_VAR} "${RNANUE_SANITIZED_FLAGS}" CACHE STRING "" FORCE)
            message(STATUS
                "Removed -stdlib=libc++ from ${RNANUE_FLAG_VAR} for GCC/libstdc++ compatibility")
        endif()
    endforeach()
endfunction()

function(rnanue_find_matching_c_compiler RNANUE_CXX_COMPILER RNANUE_OUT)
    set(RNANUE_C_COMPILER_NAMES "")
    set(RNANUE_C_COMPILER_HINTS "")
    set(RNANUE_CXX_PATHS "${RNANUE_CXX_COMPILER}")

    get_filename_component(RNANUE_CXX_REALPATH "${RNANUE_CXX_COMPILER}" REALPATH)
    if(RNANUE_CXX_REALPATH AND NOT RNANUE_CXX_REALPATH STREQUAL RNANUE_CXX_COMPILER)
        list(APPEND RNANUE_CXX_PATHS "${RNANUE_CXX_REALPATH}")
    endif()

    foreach(RNANUE_CXX_PATH IN LISTS RNANUE_CXX_PATHS)
        get_filename_component(RNANUE_CXX_DIR "${RNANUE_CXX_PATH}" DIRECTORY)
        get_filename_component(RNANUE_CXX_NAME "${RNANUE_CXX_PATH}" NAME)
        list(APPEND RNANUE_C_COMPILER_HINTS "${RNANUE_CXX_DIR}")

        if(RNANUE_CXX_NAME MATCHES "^g\\+\\+")
            string(REGEX REPLACE "^g\\+\\+" "gcc" RNANUE_GCC_NAME "${RNANUE_CXX_NAME}")
            list(APPEND RNANUE_C_COMPILER_NAMES "${RNANUE_GCC_NAME}")
        endif()
    endforeach()

    get_filename_component(RNANUE_CXX_NAME "${RNANUE_CXX_COMPILER}" NAME)
    if(RNANUE_CXX_NAME STREQUAL "c++")
        list(APPEND RNANUE_C_COMPILER_NAMES cc gcc)
    endif()

    if(NOT RNANUE_C_COMPILER_NAMES)
        set(${RNANUE_OUT} "" PARENT_SCOPE)
        return()
    endif()

    list(REMOVE_DUPLICATES RNANUE_C_COMPILER_NAMES)
    list(REMOVE_DUPLICATES RNANUE_C_COMPILER_HINTS)

    find_program(RNANUE_MATCHING_C_COMPILER
                 NAMES ${RNANUE_C_COMPILER_NAMES}
                 HINTS ${RNANUE_C_COMPILER_HINTS}
                 NO_DEFAULT_PATH
                 NO_CACHE)

    if(NOT RNANUE_MATCHING_C_COMPILER)
        find_program(RNANUE_MATCHING_C_COMPILER
                     NAMES ${RNANUE_C_COMPILER_NAMES}
                     NO_CACHE)
    endif()

    set(${RNANUE_OUT} "${RNANUE_MATCHING_C_COMPILER}" PARENT_SCOPE)
endfunction()

function(rnanue_configure_c_compiler_from_cxx RNANUE_SELECTED_CXX)
    if(DEFINED CMAKE_C_COMPILER AND NOT "${CMAKE_C_COMPILER}" STREQUAL "")
        return()
    endif()

    rnanue_find_matching_c_compiler("${RNANUE_SELECTED_CXX}" RNANUE_SELECTED_C)
    if(RNANUE_SELECTED_C)
        set(CMAKE_C_COMPILER "${RNANUE_SELECTED_C}" CACHE FILEPATH "C compiler" FORCE)
        message(STATUS "RNAnue selected matching C compiler: ${RNANUE_SELECTED_C}")
    elseif(CMAKE_HOST_APPLE)
        message(FATAL_ERROR
            "RNAnue found compatible C++ compiler '${RNANUE_SELECTED_CXX}', but could not "
            "find the matching GCC C compiler. Install the matching gcc package or configure "
            "with -DCMAKE_C_COMPILER=/path/to/gcc-${RNANUE_MINIMUM_GNU_COMPILER_VERSION}.")
    endif()
endfunction()

function(rnanue_accept_cxx_compiler RNANUE_SELECTED_CXX RNANUE_VERSION)
    set(CMAKE_CXX_COMPILER "${RNANUE_SELECTED_CXX}" CACHE FILEPATH "C++ compiler" FORCE)
    rnanue_sanitize_gnu_compiler_flags()
    rnanue_configure_c_compiler_from_cxx("${RNANUE_SELECTED_CXX}")
    message(STATUS
        "RNAnue selected C++ compiler: ${RNANUE_SELECTED_CXX} (GNU ${RNANUE_VERSION})")
endfunction()

function(rnanue_select_compatible_compilers)
    if(DEFINED CMAKE_CXX_COMPILER AND NOT "${CMAKE_CXX_COMPILER}" STREQUAL "")
        set(RNANUE_REQUESTED_CXX "${CMAKE_CXX_COMPILER}")
    else()
        set(RNANUE_REQUESTED_CXX "")
    endif()

    if(RNANUE_REQUESTED_CXX)
        rnanue_probe_cxx_compiler("${RNANUE_REQUESTED_CXX}"
                                  RNANUE_CXX_OK
                                  RNANUE_RESOLVED_CXX
                                  RNANUE_CXX_VERSION
                                  RNANUE_CXX_REASON)
        if(NOT RNANUE_CXX_OK)
            message(FATAL_ERROR
                "RNAnue cannot use the configured CMAKE_CXX_COMPILER value "
                "'${RNANUE_REQUESTED_CXX}': ${RNANUE_CXX_REASON}. "
                "Use GCC/G++ ${RNANUE_MINIMUM_GNU_COMPILER_VERSION} or newer and configure "
                "a fresh build directory with "
                "-DCMAKE_CXX_COMPILER=/path/to/g++-${RNANUE_MINIMUM_GNU_COMPILER_VERSION} "
                "-DCMAKE_C_COMPILER=/path/to/gcc-${RNANUE_MINIMUM_GNU_COMPILER_VERSION}.")
        endif()

        rnanue_accept_cxx_compiler("${RNANUE_RESOLVED_CXX}" "${RNANUE_CXX_VERSION}")
        return()
    endif()

    set(RNANUE_CXX_CANDIDATES "")
    set(RNANUE_ENV_CXX_COMPILER "")
    if(DEFINED ENV{CXX} AND NOT "$ENV{CXX}" STREQUAL "")
        set(RNANUE_ENV_CXX_COMPILER "$ENV{CXX}")
        list(APPEND RNANUE_CXX_CANDIDATES "${RNANUE_ENV_CXX_COMPILER}")
    endif()

    find_program(RNANUE_DEFAULT_CXX_COMPILER NAMES c++ g++ NO_CACHE)
    if(RNANUE_DEFAULT_CXX_COMPILER)
        list(APPEND RNANUE_CXX_CANDIDATES "${RNANUE_DEFAULT_CXX_COMPILER}")
    endif()

    if(CMAKE_HOST_APPLE)
        list(APPEND RNANUE_CXX_CANDIDATES
             /opt/homebrew/bin/g++-16
             /opt/homebrew/bin/g++-15
             /opt/homebrew/bin/g++-14
             /usr/local/bin/g++-16
             /usr/local/bin/g++-15
             /usr/local/bin/g++-14
             /opt/local/bin/g++-16
             /opt/local/bin/g++-15
             /opt/local/bin/g++-14)
    endif()

    list(APPEND RNANUE_CXX_CANDIDATES g++-16 g++-15 g++-14)
    list(REMOVE_DUPLICATES RNANUE_CXX_CANDIDATES)

    set(RNANUE_TRIED_CXX_COMPILERS "")
    set(RNANUE_ENV_CXX_REASON "")
    set(RNANUE_DEFAULT_CXX_REASON "")

    foreach(RNANUE_CXX_CANDIDATE IN LISTS RNANUE_CXX_CANDIDATES)
        rnanue_probe_cxx_compiler("${RNANUE_CXX_CANDIDATE}"
                                  RNANUE_CXX_OK
                                  RNANUE_RESOLVED_CXX
                                  RNANUE_CXX_VERSION
                                  RNANUE_CXX_REASON)

        if(RNANUE_RESOLVED_CXX)
            set(RNANUE_TRIED_CXX_NAME "${RNANUE_RESOLVED_CXX}")
        else()
            set(RNANUE_TRIED_CXX_NAME "${RNANUE_CXX_CANDIDATE}")
        endif()
        list(APPEND RNANUE_TRIED_CXX_COMPILERS
             "${RNANUE_TRIED_CXX_NAME}: ${RNANUE_CXX_REASON}")

        if(RNANUE_ENV_CXX_COMPILER
           AND RNANUE_CXX_CANDIDATE STREQUAL RNANUE_ENV_CXX_COMPILER)
            set(RNANUE_ENV_CXX_REASON "${RNANUE_CXX_REASON}")
        endif()

        if(RNANUE_DEFAULT_CXX_COMPILER
           AND RNANUE_CXX_CANDIDATE STREQUAL RNANUE_DEFAULT_CXX_COMPILER)
            set(RNANUE_DEFAULT_CXX_REASON "${RNANUE_CXX_REASON}")
        endif()

        if(RNANUE_CXX_OK)
            if(RNANUE_ENV_CXX_COMPILER
               AND NOT RNANUE_CXX_CANDIDATE STREQUAL RNANUE_ENV_CXX_COMPILER)
                message(STATUS
                    "RNAnue environment CXX is not compatible: ${RNANUE_ENV_CXX_REASON}")
            endif()
            if(RNANUE_DEFAULT_CXX_COMPILER
               AND NOT RNANUE_CXX_CANDIDATE STREQUAL RNANUE_DEFAULT_CXX_COMPILER)
                message(STATUS
                    "RNAnue default C++ compiler is not compatible: ${RNANUE_DEFAULT_CXX_REASON}")
            endif()
            rnanue_accept_cxx_compiler("${RNANUE_RESOLVED_CXX}" "${RNANUE_CXX_VERSION}")
            return()
        endif()
    endforeach()

    string(REPLACE ";" "\n  " RNANUE_TRIED_CXX_COMPILERS_TEXT
           "${RNANUE_TRIED_CXX_COMPILERS}")
    message(FATAL_ERROR
        "RNAnue requires GCC/G++ ${RNANUE_MINIMUM_GNU_COMPILER_VERSION} or newer with "
        "the required C++23 standard-library features. No compatible compiler was found.\n"
        "Tried:\n  ${RNANUE_TRIED_CXX_COMPILERS_TEXT}\n"
        "Install GCC/G++ ${RNANUE_MINIMUM_GNU_COMPILER_VERSION} or newer, put it on PATH, "
        "or configure a fresh build directory with "
        "-DCMAKE_CXX_COMPILER=/path/to/g++-${RNANUE_MINIMUM_GNU_COMPILER_VERSION} "
        "-DCMAKE_C_COMPILER=/path/to/gcc-${RNANUE_MINIMUM_GNU_COMPILER_VERSION}.")
endfunction()
