set(RNANUE_DEPENDENCY_PROVIDER "AUTO" CACHE STRING
    "Default dependency provider: AUTO tries installed packages first, SYSTEM requires installed packages, BUNDLED builds pinned bundled fallbacks")
set_property(CACHE RNANUE_DEPENDENCY_PROVIDER PROPERTY STRINGS AUTO SYSTEM BUNDLED)

set(RNANUE_DEPENDENCY_PREFIX "" CACHE PATH
    "Preferred installed dependency prefix. When unset, Conda-style prefixes from CMAKE_PREFIX_PATH or CONDA_PREFIX are detected automatically.")
option(RNANUE_PREFER_ACTIVE_PREFIX_LIBS
    "Prefer common runtime libraries from an active Conda-style dependency prefix to avoid mixed unsafe RPATHs"
    ON)

function(rnanue_normalize_provider VAR DESCRIPTION)
    if(NOT DEFINED ${VAR} OR "${${VAR}}" STREQUAL "")
        set(${VAR} "${RNANUE_DEPENDENCY_PROVIDER}" CACHE STRING "${DESCRIPTION}" FORCE)
    endif()

    string(TOUPPER "${${VAR}}" RNANUE_PROVIDER_VALUE)
    set(${VAR} "${RNANUE_PROVIDER_VALUE}" CACHE STRING "${DESCRIPTION}" FORCE)
    set_property(CACHE ${VAR} PROPERTY STRINGS AUTO SYSTEM BUNDLED)

    if(NOT RNANUE_PROVIDER_VALUE IN_LIST RNANUE_PROVIDER_VALUES)
        message(FATAL_ERROR
            "Invalid ${VAR}='${RNANUE_PROVIDER_VALUE}'. Use AUTO, SYSTEM, or BUNDLED.")
    endif()

    set(${VAR} "${RNANUE_PROVIDER_VALUE}" PARENT_SCOPE)
endfunction()

set(RNANUE_PROVIDER_VALUES AUTO SYSTEM BUNDLED)
rnanue_normalize_provider(RNANUE_DEPENDENCY_PROVIDER
    "Default dependency provider: AUTO tries installed packages first, SYSTEM requires installed packages, BUNDLED builds pinned bundled fallbacks")

function(rnanue_configure_dependency_provider VAR DESCRIPTION)
    rnanue_normalize_provider(${VAR} "${DESCRIPTION}")
    set(${VAR} "${${VAR}}" PARENT_SCOPE)
endfunction()

function(rnanue_provider_allows_system VAR OUT_VAR)
    if("${${VAR}}" STREQUAL "AUTO" OR "${${VAR}}" STREQUAL "SYSTEM")
        set(${OUT_VAR} TRUE PARENT_SCOPE)
    else()
        set(${OUT_VAR} FALSE PARENT_SCOPE)
    endif()
endfunction()

function(rnanue_provider_allows_bundled VAR OUT_VAR)
    if("${${VAR}}" STREQUAL "AUTO" OR "${${VAR}}" STREQUAL "BUNDLED")
        set(${OUT_VAR} TRUE PARENT_SCOPE)
    else()
        set(${OUT_VAR} FALSE PARENT_SCOPE)
    endif()
endfunction()

function(rnanue_path_is_conda_like PATH_VALUE OUT_VAR)
    if(PATH_VALUE MATCHES "/(conda|miniconda|miniconda3|miniforge|miniforge3|anaconda|anaconda3|mambaforge|envs)(/|$)")
        set(${OUT_VAR} TRUE PARENT_SCOPE)
    else()
        set(${OUT_VAR} FALSE PARENT_SCOPE)
    endif()
endfunction()

function(rnanue_detect_active_dependency_prefix OUT_VAR)
    if(RNANUE_DEPENDENCY_PREFIX)
        set(${OUT_VAR} "${RNANUE_DEPENDENCY_PREFIX}" PARENT_SCOPE)
        return()
    endif()

    set(RNANUE_PREFIX_CANDIDATES ${CMAKE_PREFIX_PATH})
    if(DEFINED ENV{CONDA_PREFIX} AND NOT "$ENV{CONDA_PREFIX}" STREQUAL "")
        list(APPEND RNANUE_PREFIX_CANDIDATES "$ENV{CONDA_PREFIX}")
    endif()

    foreach(RNANUE_PREFIX IN LISTS RNANUE_PREFIX_CANDIDATES)
        if(NOT RNANUE_PREFIX)
            continue()
        endif()

        rnanue_path_is_conda_like("${RNANUE_PREFIX}" RNANUE_IS_CONDA_PREFIX)
        if(RNANUE_IS_CONDA_PREFIX AND EXISTS "${RNANUE_PREFIX}/lib")
            set(${OUT_VAR} "${RNANUE_PREFIX}" PARENT_SCOPE)
            return()
        endif()
    endforeach()

    set(${OUT_VAR} "" PARENT_SCOPE)
endfunction()

function(rnanue_find_prefix_library OUT_VAR PREFIX)
    find_library(RNANUE_PREFIX_LIBRARY
        NAMES ${ARGN}
        PATHS "${PREFIX}/lib" "${PREFIX}/lib64"
        NO_DEFAULT_PATH
        NO_CACHE
    )
    set(${OUT_VAR} "${RNANUE_PREFIX_LIBRARY}" PARENT_SCOPE)
endfunction()

function(rnanue_find_prefix_include OUT_VAR PREFIX)
    find_path(RNANUE_PREFIX_INCLUDE
        NAMES ${ARGN}
        PATHS "${PREFIX}/include" "${PREFIX}/include/libpng16"
        NO_DEFAULT_PATH
        NO_CACHE
    )
    set(${OUT_VAR} "${RNANUE_PREFIX_INCLUDE}" PARENT_SCOPE)
endfunction()

function(rnanue_set_prefix_cache_path VAR VALUE DOC)
    if(NOT VALUE)
        return()
    endif()

    set(RNANUE_SHOULD_SET FALSE)
    if(RNANUE_DEPENDENCY_PREFIX)
        set(RNANUE_SHOULD_SET TRUE)
    elseif(NOT DEFINED ${VAR} OR "${${VAR}}" STREQUAL "" OR "${${VAR}}" MATCHES "^/usr/" OR "${${VAR}}" MATCHES "^/lib/")
        set(RNANUE_SHOULD_SET TRUE)
    endif()

    if(RNANUE_SHOULD_SET)
        set(${VAR} "${VALUE}" CACHE PATH "${DOC}" FORCE)
    endif()
endfunction()

function(rnanue_prefer_active_prefix_libraries PREFIX)
    if(NOT RNANUE_PREFER_ACTIVE_PREFIX_LIBS OR NOT PREFIX)
        return()
    endif()

    message(STATUS "RNAnue prefers common runtime libraries from dependency prefix: ${PREFIX}")

    list(PREPEND CMAKE_PREFIX_PATH "${PREFIX}")
    list(PREPEND CMAKE_LIBRARY_PATH "${PREFIX}/lib" "${PREFIX}/lib64")
    list(PREPEND CMAKE_INCLUDE_PATH "${PREFIX}/include")
    set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH}" PARENT_SCOPE)
    set(CMAKE_LIBRARY_PATH "${CMAKE_LIBRARY_PATH}" PARENT_SCOPE)
    set(CMAKE_INCLUDE_PATH "${CMAKE_INCLUDE_PATH}" PARENT_SCOPE)

    set(ZLIB_ROOT "${PREFIX}" CACHE PATH "Preferred zlib prefix" FORCE)
    set(BZip2_ROOT "${PREFIX}" CACHE PATH "Preferred BZip2 prefix" FORCE)
    set(PNG_ROOT "${PREFIX}" CACHE PATH "Preferred libpng prefix" FORCE)
    set(TBB_ROOT "${PREFIX}" CACHE PATH "Preferred oneTBB prefix" FORCE)

    rnanue_find_prefix_include(RNANUE_PREFIX_ZLIB_INCLUDE "${PREFIX}" zlib.h)
    rnanue_find_prefix_library(RNANUE_PREFIX_ZLIB_LIBRARY "${PREFIX}" z zlib libz.so.1)
    rnanue_set_prefix_cache_path(ZLIB_INCLUDE_DIR "${RNANUE_PREFIX_ZLIB_INCLUDE}" "Preferred zlib include directory")
    rnanue_set_prefix_cache_path(ZLIB_LIBRARY "${RNANUE_PREFIX_ZLIB_LIBRARY}" "Preferred zlib library")

    rnanue_find_prefix_include(RNANUE_PREFIX_BZIP2_INCLUDE "${PREFIX}" bzlib.h)
    rnanue_find_prefix_library(RNANUE_PREFIX_BZIP2_LIBRARY "${PREFIX}" bz2 bzip2 libbz2.so.1.0)
    rnanue_set_prefix_cache_path(BZIP2_INCLUDE_DIR "${RNANUE_PREFIX_BZIP2_INCLUDE}" "Preferred BZip2 include directory")
    rnanue_set_prefix_cache_path(BZIP2_LIBRARY "${RNANUE_PREFIX_BZIP2_LIBRARY}" "Preferred BZip2 library")
    rnanue_set_prefix_cache_path(BZIP2_LIBRARY_RELEASE "${RNANUE_PREFIX_BZIP2_LIBRARY}" "Preferred BZip2 release library")

    rnanue_find_prefix_include(RNANUE_PREFIX_PNG_INCLUDE "${PREFIX}" png.h)
    rnanue_find_prefix_library(RNANUE_PREFIX_PNG_LIBRARY "${PREFIX}" png16 png libpng16.so.16)
    rnanue_set_prefix_cache_path(PNG_PNG_INCLUDE_DIR "${RNANUE_PREFIX_PNG_INCLUDE}" "Preferred libpng include directory")
    rnanue_set_prefix_cache_path(PNG_LIBRARY "${RNANUE_PREFIX_PNG_LIBRARY}" "Preferred libpng library")

    rnanue_find_prefix_include(RNANUE_PREFIX_TBB_INCLUDE "${PREFIX}" tbb/tbb.h oneapi/tbb.h)
    rnanue_find_prefix_library(RNANUE_PREFIX_TBB_LIBRARY "${PREFIX}" tbb libtbb.so.12)
    set(RNANUE_PREFIX_TBB_INCLUDE_DIR "${RNANUE_PREFIX_TBB_INCLUDE}" CACHE PATH "Preferred oneTBB include directory" FORCE)
    set(RNANUE_PREFIX_TBB_LIBRARY "${RNANUE_PREFIX_TBB_LIBRARY}" CACHE FILEPATH "Preferred oneTBB library" FORCE)
endfunction()

function(rnanue_try_prefix_tbb OUT_VAR)
    set(${OUT_VAR} FALSE PARENT_SCOPE)

    if(TARGET TBB::tbb OR NOT RNANUE_PREFIX_TBB_LIBRARY)
        return()
    endif()

    add_library(TBB::tbb SHARED IMPORTED GLOBAL)
    set_target_properties(TBB::tbb PROPERTIES
        IMPORTED_LOCATION "${RNANUE_PREFIX_TBB_LIBRARY}"
    )

    if(RNANUE_PREFIX_TBB_INCLUDE_DIR)
        set_target_properties(TBB::tbb PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${RNANUE_PREFIX_TBB_INCLUDE_DIR}"
        )
    endif()

    set(TBB_FOUND TRUE PARENT_SCOPE)
    set(${OUT_VAR} TRUE PARENT_SCOPE)
endfunction()

function(rnanue_classify_dependency_origin OUT_VAR)
    set(RNANUE_PATHS ${ARGN})
    set(RNANUE_ORIGIN "unknown")

    foreach(RNANUE_PATH IN LISTS RNANUE_PATHS)
        if(NOT RNANUE_PATH)
            continue()
        endif()

        if(RNANUE_PATH MATCHES "^${CMAKE_BINARY_DIR}/")
            set(RNANUE_ORIGIN "bundled")
            break()
        elseif(RNANUE_PATH MATCHES "/(conda|miniconda|miniconda3|miniforge|miniforge3|anaconda|anaconda3|mambaforge|envs)/")
            set(RNANUE_ORIGIN "conda")
            break()
        elseif(RNANUE_PATH MATCHES "^/opt/homebrew/" OR RNANUE_PATH MATCHES "^/usr/local/")
            set(RNANUE_ORIGIN "homebrew")
        elseif(RNANUE_PATH MATCHES "^/opt/local/")
            set(RNANUE_ORIGIN "macports")
        elseif(RNANUE_PATH MATCHES "^/usr/")
            if(RNANUE_ORIGIN STREQUAL "unknown")
                set(RNANUE_ORIGIN "system")
            endif()
        elseif(RNANUE_ORIGIN STREQUAL "unknown")
            set(RNANUE_ORIGIN "custom")
        endif()
    endforeach()

    set(${OUT_VAR} "${RNANUE_ORIGIN}" PARENT_SCOPE)
endfunction()

function(rnanue_warn_if_shadowing_implicit NAME)
    if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
        return()
    endif()

    foreach(RNANUE_PATH IN LISTS ARGN)
        if(NOT RNANUE_PATH OR RNANUE_PATH MATCHES "^/usr/lib" OR RNANUE_PATH MATCHES "^/lib")
            continue()
        endif()

        if(IS_DIRECTORY "${RNANUE_PATH}")
            set(RNANUE_LIB_DIR "${RNANUE_PATH}")
        else()
            get_filename_component(RNANUE_LIB_DIR "${RNANUE_PATH}" DIRECTORY)
        endif()

        if(EXISTS "${RNANUE_LIB_DIR}/libz.so.1")
            message(WARNING
                "${NAME} uses non-system library directory '${RNANUE_LIB_DIR}', which contains "
                "libz.so.1 and may shadow the implicit system zlib. Prefer a consistent "
                "dependency prefix via CMAKE_PREFIX_PATH/ZLIB_ROOT, or use "
                "-DRNANUE_DEPENDENCY_PROVIDER=SYSTEM for no-download system-only builds.")
        endif()
    endforeach()
endfunction()

function(rnanue_record_dependency NAME PROVIDER INCLUDE_PATHS LIBRARIES)
    rnanue_classify_dependency_origin(RNANUE_ORIGIN ${INCLUDE_PATHS} ${LIBRARIES})
    rnanue_warn_if_shadowing_implicit("${NAME}" ${LIBRARIES})

    string(REPLACE ";" ", " RNANUE_INCLUDE_TEXT "${INCLUDE_PATHS}")
    string(REPLACE ";" ", " RNANUE_LIBRARY_TEXT "${LIBRARIES}")
    set(RNANUE_LINE
        "  ${NAME}: provider=${PROVIDER} | origin=${RNANUE_ORIGIN} | includes=${RNANUE_INCLUDE_TEXT} | libs=${RNANUE_LIBRARY_TEXT}")
    set_property(GLOBAL APPEND PROPERTY RNANUE_DEPENDENCY_SUMMARY "${RNANUE_LINE}")
endfunction()

function(rnanue_print_dependency_summary)
    get_property(RNANUE_SUMMARY GLOBAL PROPERTY RNANUE_DEPENDENCY_SUMMARY)
    if(RNANUE_SUMMARY)
        message(STATUS "RNAnue dependency summary:")
        foreach(RNANUE_LINE IN LISTS RNANUE_SUMMARY)
            message(STATUS "${RNANUE_LINE}")
        endforeach()
    endif()
endfunction()
