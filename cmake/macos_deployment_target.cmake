# Older GCC releases infer macOS 16 from Darwin 25. Apple's assembler instead
# uses macOS 26, producing conflicting deployment versions. Give native builds
# an explicit default on macOS 26+; retain user/toolchain targets and leave
# cross-compilation and older hosts unchanged. Include after project() so
# CMake has already loaded toolchain and MACOSX_DEPLOYMENT_TARGET settings.
if(APPLE AND NOT CMAKE_CROSSCOMPILING AND NOT CMAKE_OSX_DEPLOYMENT_TARGET)
    execute_process(
        COMMAND /usr/bin/sw_vers -productVersion
        RESULT_VARIABLE RNANUE_MACOS_VERSION_RESULT
        OUTPUT_VARIABLE RNANUE_MACOS_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(
        RNANUE_MACOS_VERSION_RESULT EQUAL 0
        AND RNANUE_MACOS_VERSION MATCHES "^([0-9]+)\\."
        AND RNANUE_MACOS_VERSION VERSION_GREATER_EQUAL "26.0"
    )
        set(CMAKE_OSX_DEPLOYMENT_TARGET
            "${CMAKE_MATCH_1}.0"
            CACHE STRING
            "Minimum macOS version to target"
            FORCE
        )
        message(
            STATUS
            "RNAnue default macOS deployment target: ${CMAKE_OSX_DEPLOYMENT_TARGET}"
        )
    endif()
endif()
