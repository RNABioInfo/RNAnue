cmake_minimum_required(VERSION 3.24)
project(RNAnueStarPackaging VERSION 1.0.0 LANGUAGES C CXX)
include(CTest)
include("${CMAKE_CURRENT_LIST_DIR}/../../cmake/star.cmake")
add_executable(star_locator_probe locator.cpp)
rnanue_add_star(star_locator_probe)
install(TARGETS star_locator_probe RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(NAME StarRuntimeLicense COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_LIST_DIR}/runtime_license.py" --cmake "${CMAKE_COMMAND}")
if(RNANUE_BUILD_STAR)
    set(static_argument)
    if(RNANUE_STAR_STATIC)
        set(static_argument --static)
    endif()
    add_test(NAME StarSmoke COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/smoke.py"
        --star "$<TARGET_FILE:STAR::STAR>" --locator "$<TARGET_FILE:star_locator_probe>")
    add_test(NAME StarPackage COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_LIST_DIR}/package.py"
        --build "${CMAKE_BINARY_DIR}" --cmake "${CMAKE_COMMAND}"
        --libexec "${CMAKE_INSTALL_LIBEXECDIR}" --bindir "${CMAKE_INSTALL_BINDIR}"
        --datadir "${CMAKE_INSTALL_DATADIR}" ${static_argument})
    set_tests_properties(StarSmoke StarPackage PROPERTIES TIMEOUT 120)
else()
    add_test(NAME StarDisabled COMMAND star_locator_probe --expect-disabled)
endif()
set(CPACK_PACKAGE_NAME RNAnueStarPackaging)
set(CPACK_PACKAGE_CONTACT "RNAnue maintainers")
set(CPACK_GENERATOR TGZ)
include(CPack)
