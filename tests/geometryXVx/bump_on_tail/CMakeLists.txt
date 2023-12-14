# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(NAME TestSimulationBumpontailFFT
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_bumpontail.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:bumpontail_fft>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "fft")
set_property(TEST TestSimulationBumpontailFFT PROPERTY TIMEOUT 200)
set_property(TEST TestSimulationBumpontailFFT PROPERTY COST 100)

add_test(NAME TestSimulationBumpontailRestartFFT_XVx
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_bumpontail_restart.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:bumpontail_fft>"
        "restart")
set_property(TEST TestSimulationBumpontailRestartFFT_XVx PROPERTY TIMEOUT 200)
set_property(TEST TestSimulationBumpontailRestartFFT_XVx PROPERTY COST 100)
set_tests_properties(TestSimulationBumpontailRestartFFT_XVx PROPERTIES LABELS "ResultsNotTested" )

add_test(NAME TestSimulationBumpontailFemUniform_xperiod_vx
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_bumpontail.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:bumpontail_fem_uniform_xperiod_vx>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "fem")
set_property(TEST TestSimulationBumpontailFemUniform_xperiod_vx PROPERTY TIMEOUT 200)
set_property(TEST TestSimulationBumpontailFemUniform_xperiod_vx PROPERTY COST 100)

add_test(NAME TestSimulationBumpontailRestartFemUniform_xperiod_vx
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_bumpontail_restart.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:bumpontail_fem_uniform_xperiod_vx>"
        "restart")
set_property(TEST TestSimulationBumpontailRestartFemUniform_xperiod_vx PROPERTY TIMEOUT 200)
set_property(TEST TestSimulationBumpontailRestartFemUniform_xperiod_vx PROPERTY COST 100)
set_tests_properties(TestSimulationBumpontailRestartFemUniform_xperiod_vx PROPERTIES LABELS "ResultsNotTested" )
