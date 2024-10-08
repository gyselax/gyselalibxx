# SPDX-License-Identifier: MIT

find_package(Python3 REQUIRED COMPONENTS Interpreter)

if (${ACTIVATE_RESTART_TESTS})
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

    add_test(NAME TestSimulationSheathRestart_xperiod_vx_user_defined_grid
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_sheath_restart_user_defined_grid.sh"
            "${PROJECT_SOURCE_DIR}"
            "$<TARGET_FILE:sheath_xperiod_vx_user_defined_grid>"
            "$<TARGET_FILE:Python3::Interpreter>"
            "restart"
            "${RELATIVE_RESTART_TOLERANCE}"
            "${ABSOLUTE_RESTART_TOLERANCE}")
    set_property(TEST TestSimulationSheathRestart_xperiod_vx_user_defined_grid PROPERTY TIMEOUT 200)
    set_property(TEST TestSimulationSheathRestart_xperiod_vx_user_defined_grid PROPERTY COST 100)
endif()
