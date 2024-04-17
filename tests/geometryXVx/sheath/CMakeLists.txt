# SPDX-License-Identifier: MIT

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(NAME TestSimulationSheathRestart_xperiod_vx
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_sheath_restart.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:sheath_xperiod_vx>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "restart"
        "${RELATIVE_RESTART_TOLERANCE}"
        "${ABSOLUTE_RESTART_TOLERANCE}")
set_property(TEST TestSimulationSheathRestart_xperiod_vx PROPERTY TIMEOUT 200)
set_property(TEST TestSimulationSheathRestart_xperiod_vx PROPERTY COST 100)
set_tests_properties(TestSimulationSheathRestart_xperiod_vx PROPERTIES LABEL "ResultsNotTested" )
