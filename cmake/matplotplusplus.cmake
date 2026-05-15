option(CPM_USE_LOCAL_PACKAGES
       "Try `find_package` before downloading dependencies" ON)

set(MATPLOTPP_BUILD_EXAMPLES OFF CACHE BOOL "Build Matplot++ examples" FORCE)
set(MATPLOTPP_BUILD_TESTS OFF CACHE BOOL "Build Matplot++ tests" FORCE)
set(MATPLOTPP_BUILD_INSTALLER OFF CACHE BOOL "Build Matplot++ installer target" FORCE)

CPMAddPackage(
    NAME matplotplusplus
    GITHUB_REPOSITORY alandefreitas/matplotplusplus
    GIT_TAG origin/master # or whatever tag you want
)

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
