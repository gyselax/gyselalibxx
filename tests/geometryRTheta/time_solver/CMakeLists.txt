# SPDX-License-Identifier: MIT

# foreach(PREDCORR "RK2_PREDCORR" "EXPLICIT_PREDCORR" "IMPLICIT_PREDCORR")
#     set(test_name "bsl_${PREDCORR}_test")
#     add_executable("${test_name}"
#         ../../main.cpp
#         bsl_predcorr_test.cpp
#     )
#     target_link_libraries("${test_name}"
#         PUBLIC
#         GTest::gtest
#         GTest::gmock
#         DDC::core
#         DDC::pdi

#         paraconf::paraconf

#         PDI::pdi

#         gslx::advection
#         gslx::diocotron_init_eq
#         gslx::geometry_RTheta
#         gslx::interpolation
#         gslx::io
#         gslx::paraconfpp
#         gslx::poisson_RTheta
#         gslx::quadrature
#         gslx::simulation_utils
#         gslx::time_integration_RTheta
#         gslx::utils
#     )
#     target_compile_definitions("${test_name}" PUBLIC -D${PREDCORR} )
#     gtest_discover_tests("${test_name}" DISCOVERY_MODE PRE_TEST)
# endforeach()


find_package(Python3 REQUIRED COMPONENTS Interpreter)

function (diocotron_executable PREDCORR_METHOD TIME_METHOD )
    set(test_name "TestSimulationDiocotron_${PREDCORR_METHOD}_${TIME_METHOD}")

    add_test(NAME "${test_name}"
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_diocotron.sh"
            "${PROJECT_SOURCE_DIR}"
            "$<TARGET_FILE:diocotron_${PREDCORR_METHOD}_${TIME_METHOD}>"
            "$<TARGET_FILE:Python3::Interpreter>")
            # "fft")
            
    # set_property(TEST "${test_name}" PROPERTY TIMEOUT 2000)
    # set_property(TEST "${test_name}" PROPERTY COST 100)
endfunction()

# Select the parameters of the simulation:
diocotron_executable(PREDCORR RK3_METHOD)
diocotron_executable(EXPLICIT_PREDCORR EULER_METHOD)
diocotron_executable(IMPLICIT_PREDCORR EULER_METHOD)
