cmake_minimum_required(VERSION 3.24)

foreach(RNANUE_REQUIRED_VAR
        RNANUE_BOOST_STEP
        RNANUE_BOOST_SOURCE_DIR
        RNANUE_BOOST_INSTALL_DIR
        RNANUE_BOOST_CXX_COMPILER
        RNANUE_BOOST_TOOLSET)
    if(NOT DEFINED ${RNANUE_REQUIRED_VAR} OR "${${RNANUE_REQUIRED_VAR}}" STREQUAL "")
        message(FATAL_ERROR "${RNANUE_REQUIRED_VAR} is required")
    endif()
endforeach()

if(NOT DEFINED RNANUE_BOOST_CXXFLAGS)
    set(RNANUE_BOOST_CXXFLAGS "")
endif()

function(rnanue_run_boost_step)
    cmake_parse_arguments(RNANUE_STEP "" "WORKING_DIRECTORY" "COMMAND" ${ARGN})

    execute_process(
        COMMAND ${RNANUE_STEP_COMMAND}
        WORKING_DIRECTORY "${RNANUE_STEP_WORKING_DIRECTORY}"
        RESULT_VARIABLE RNANUE_STEP_RESULT
    )

    if(NOT RNANUE_STEP_RESULT EQUAL 0)
        message(FATAL_ERROR "Boost ${RNANUE_BOOST_STEP} step failed with exit code ${RNANUE_STEP_RESULT}")
    endif()
endfunction()

if(RNANUE_BOOST_STEP STREQUAL "bootstrap")
    set(RNANUE_B2_ENGINE "${RNANUE_BOOST_SOURCE_DIR}/tools/build/src/engine/b2")

    set(RNANUE_B2_ENGINE_COMMAND
        "${CMAKE_COMMAND}" -E env
        "CXX=${RNANUE_BOOST_CXX_COMPILER}"
        ./tools/build/src/engine/build.sh
        "--cxx=${RNANUE_BOOST_CXX_COMPILER}"
        "--cxxflags=${RNANUE_BOOST_CXXFLAGS}"
        "${RNANUE_BOOST_TOOLSET}"
    )

    rnanue_run_boost_step(
        COMMAND ${RNANUE_B2_ENGINE_COMMAND}
        WORKING_DIRECTORY "${RNANUE_BOOST_SOURCE_DIR}"
    )

    rnanue_run_boost_step(
        COMMAND "${CMAKE_COMMAND}" -E copy "${RNANUE_B2_ENGINE}" "${RNANUE_BOOST_SOURCE_DIR}/b2"
        WORKING_DIRECTORY "${RNANUE_BOOST_SOURCE_DIR}"
    )

    rnanue_run_boost_step(
        COMMAND "${CMAKE_COMMAND}" -E env
            "CXX=${RNANUE_BOOST_CXX_COMPILER}"
            ./bootstrap.sh
            "--prefix=${RNANUE_BOOST_INSTALL_DIR}"
            "--with-toolset=${RNANUE_BOOST_TOOLSET}"
            --with-libraries=program_options
            "--with-bjam=${RNANUE_BOOST_SOURCE_DIR}/b2"
        WORKING_DIRECTORY "${RNANUE_BOOST_SOURCE_DIR}"
    )
elseif(RNANUE_BOOST_STEP STREQUAL "build")
    set(RNANUE_B2_COMMAND
        "${CMAKE_COMMAND}" -E env
        "CXX=${RNANUE_BOOST_CXX_COMPILER}"
        ./b2 install
        "--cxx=${RNANUE_BOOST_CXX_COMPILER}"
        "toolset=${RNANUE_BOOST_TOOLSET}"
        link=static
        variant=release
        threading=multi
        runtime-link=static
    )

    if(NOT "${RNANUE_BOOST_CXXFLAGS}" STREQUAL "")
        list(APPEND RNANUE_B2_COMMAND "cxxflags=${RNANUE_BOOST_CXXFLAGS}")
    endif()

    rnanue_run_boost_step(
        COMMAND ${RNANUE_B2_COMMAND}
        WORKING_DIRECTORY "${RNANUE_BOOST_SOURCE_DIR}"
    )
else()
    message(FATAL_ERROR "Unknown RNANUE_BOOST_STEP='${RNANUE_BOOST_STEP}'")
endif()
