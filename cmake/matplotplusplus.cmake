include(${CMAKE_CURRENT_LIST_DIR}/dependency_providers.cmake)

option(CPM_USE_LOCAL_PACKAGES
       "Try `find_package` before downloading dependencies" ON)

rnanue_configure_dependency_provider(RNANUE_MATPLOT_PROVIDER
    "Matplot++ provider: AUTO tries installed Matplot++ first, SYSTEM requires it, BUNDLED downloads a pinned Matplot++ fallback")
rnanue_provider_allows_system(RNANUE_MATPLOT_PROVIDER RNANUE_MATPLOT_ALLOW_SYSTEM)
rnanue_provider_allows_bundled(RNANUE_MATPLOT_PROVIDER RNANUE_MATPLOT_ALLOW_BUNDLED)

set(MATPLOTPP_BUILD_EXAMPLES OFF CACHE BOOL "Build Matplot++ examples" FORCE)
set(MATPLOTPP_BUILD_TESTS OFF CACHE BOOL "Build Matplot++ tests" FORCE)
set(MATPLOTPP_BUILD_INSTALLER OFF CACHE BOOL "Build Matplot++ installer target" FORCE)

if(RNANUE_MATPLOT_ALLOW_SYSTEM)
    find_package(matplotplusplus CONFIG QUIET)
endif()

if(NOT TARGET matplot AND NOT TARGET matplotplusplus::matplot AND NOT TARGET matplotplusplus::matplotplusplus)
    if(RNANUE_MATPLOT_PROVIDER STREQUAL "SYSTEM")
        message(FATAL_ERROR
            "RNANUE_MATPLOT_PROVIDER=SYSTEM requires installed Matplot++. Install Matplot++, "
            "set CMAKE_PREFIX_PATH to its prefix, or configure with "
            "-DRNANUE_MATPLOT_PROVIDER=AUTO or BUNDLED.")
    elseif(RNANUE_MATPLOT_ALLOW_BUNDLED)
        CPMAddPackage(
            NAME matplotplusplus
            GITHUB_REPOSITORY alandefreitas/matplotplusplus
            GIT_TAG 2ccbb7ff86295bb4d3aca2938dce09ca3accb1ef
        )
    endif()
endif()

if(NOT TARGET matplot)
    if(TARGET matplotplusplus::matplot)
        add_library(matplot ALIAS matplotplusplus::matplot)
    elseif(TARGET matplotplusplus::matplotplusplus)
        add_library(matplot ALIAS matplotplusplus::matplotplusplus)
    else()
        message(FATAL_ERROR "Matplot++ target 'matplot' was not created.")
    endif()
endif()

rnanue_record_dependency("Matplot++" "${RNANUE_MATPLOT_PROVIDER}" "" "matplot")

find_package(ZLIB QUIET)
find_package(PNG QUIET)

if(TARGET cimg)
    if(TARGET ZLIB::ZLIB)
        target_link_libraries(cimg INTERFACE ZLIB::ZLIB)
    endif()

    if(TARGET PNG::PNG)
        target_link_libraries(cimg INTERFACE PNG::PNG)
    endif()

    foreach(RNANUE_MATPLOT_INCLUDE_DIR IN LISTS ZLIB_INCLUDE_DIRS PNG_INCLUDE_DIRS)
        if(RNANUE_MATPLOT_INCLUDE_DIR AND EXISTS "${RNANUE_MATPLOT_INCLUDE_DIR}")
            target_include_directories(cimg INTERFACE "${RNANUE_MATPLOT_INCLUDE_DIR}")

            if(TARGET matplot)
                target_compile_options(matplot PRIVATE
                    "$<$<COMPILE_LANGUAGE:CXX>:-I${RNANUE_MATPLOT_INCLUDE_DIR}>")
            endif()
        endif()
    endforeach()
endif()
