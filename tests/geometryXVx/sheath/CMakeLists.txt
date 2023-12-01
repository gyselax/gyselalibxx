# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

add_test(NAME TestSimulationSheathRestart_xperiod_vx
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_sheath_restart.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:sheath_xperiod_vx>"
        "restart")
set_property(TEST TestSimulationSheathRestart_xperiod_vx PROPERTY TIMEOUT 200)
set_tests_properties(TestSimulationSheathRestart_xperiod_vx PROPERTIES LABEL "ResultsNotTested" )
