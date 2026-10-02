include_guard(GLOBAL)

# Keep toolchain notices out of the source tree. Distribution copyright files
# may contain additional notices; copy them intact rather than extracting text.
function(rnanue_stage_star_runtime_license destination)
    set(license "")
    set(RNANUE_STAR_GCC_RUNTIME_LICENSE
        ""
        CACHE FILEPATH
        "GCC COPYING.RUNTIME or toolchain copyright file containing the runtime exception"
    )
    if(RNANUE_STAR_GCC_RUNTIME_LICENSE)
        set(license "${RNANUE_STAR_GCC_RUNTIME_LICENSE}")
    else()
        get_filename_component(compiler "${CMAKE_CXX_COMPILER}" REALPATH)
        get_filename_component(compiler_bin "${compiler}" DIRECTORY)
        get_filename_component(prefix "${compiler_bin}" DIRECTORY)
        string(REGEX MATCH "^[0-9]+" major "${CMAKE_CXX_COMPILER_VERSION}")
        set(candidates
            "${prefix}/COPYING.RUNTIME"
            "${prefix}/share/doc/gcc-${major}-base/copyright"
            "${prefix}/share/doc/gcc-${major}/COPYING.RUNTIME"
            "${prefix}/share/doc/gcc-${major}/copyright"
            "${prefix}/share/doc/gcc${major}/COPYING.RUNTIME"
            "${prefix}/share/licenses/gcc/COPYING.RUNTIME"
            "${prefix}/share/licenses/gcc-libs/COPYING.RUNTIME"
            "${prefix}/share/licenses/libgcc/COPYING.RUNTIME"
            "${prefix}/share/doc/gcc/COPYING.RUNTIME"
        )
        foreach(candidate IN LISTS candidates)
            if(EXISTS "${candidate}" AND NOT IS_DIRECTORY "${candidate}")
                file(READ "${candidate}" text)
                if(text MATCHES "GCC RUNTIME LIBRARY EXCEPTION")
                    set(license "${candidate}")
                    break()
                endif()
            endif()
        endforeach()
    endif()
    if(NOT license OR NOT EXISTS "${license}" OR IS_DIRECTORY "${license}")
        message(
            FATAL_ERROR
            "Cannot find the GCC runtime exception for ${CMAKE_CXX_COMPILER}. "
            "Install the toolchain's licence/documentation files or set "
            "RNANUE_STAR_GCC_RUNTIME_LICENSE to its COPYING.RUNTIME or copyright file. "
            "The notice is staged in the build tree; do not add it to the repository."
        )
    endif()
    file(READ "${license}" text)
    if(
        NOT text MATCHES "GCC RUNTIME LIBRARY EXCEPTION"
        OR NOT text MATCHES "No Weakening of GCC Copyleft"
    )
        message(
            FATAL_ERROR
            "RNANUE_STAR_GCC_RUNTIME_LICENSE does not contain the GCC runtime exception: ${license}"
        )
    endif()
    get_filename_component(directory "${destination}" DIRECTORY)
    file(MAKE_DIRECTORY "${directory}")
    configure_file("${license}" "${destination}" COPYONLY)
    message(STATUS "STAR GCC runtime notice: ${license}")
endfunction()
