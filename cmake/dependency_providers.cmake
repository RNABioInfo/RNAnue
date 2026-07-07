set(RNANUE_DEPENDENCY_PROVIDER "AUTO" CACHE STRING
    "Default dependency provider: AUTO tries installed packages first, SYSTEM requires installed packages, BUNDLED builds pinned bundled fallbacks")
set_property(CACHE RNANUE_DEPENDENCY_PROVIDER PROPERTY STRINGS AUTO SYSTEM BUNDLED)

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
