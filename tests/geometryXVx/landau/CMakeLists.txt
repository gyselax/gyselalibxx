# SPDX-License-Identifier: MIT

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(NAME TestSimulationLandauFFT_XVx
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_landau.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:landau_fft>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "fft")
set_property(TEST TestSimulationLandauFFT_XVx PROPERTY TIMEOUT 200)
set_property(TEST TestSimulationLandauFFT_XVx PROPERTY COST 100)

if (${ACTIVATE_RESTART_TESTS})
    add_test(NAME TestSimulationLandauRestartFFT_XVx
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_landau_restart.sh"
            "${PROJECT_SOURCE_DIR}"
            "$<TARGET_FILE:landau_fft>"
            "$<TARGET_FILE:Python3::Interpreter>"
            "restart"
            "${RELATIVE_RESTART_TOLERANCE}"
            "${ABSOLUTE_RESTART_TOLERANCE}")
    set_property(TEST TestSimulationLandauRestartFFT_XVx PROPERTY TIMEOUT 200)
    set_property(TEST TestSimulationLandauRestartFFT_XVx PROPERTY COST 100)
endif()

add_test(NAME TestSimulationLandauFemUniform_xperiod_vx
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_landau.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:landau_fem_uniform_xperiod_vx>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "fem")
set_property(TEST TestSimulationLandauFemUniform_xperiod_vx PROPERTY TIMEOUT 200)
set_property(TEST TestSimulationLandauFemUniform_xperiod_vx PROPERTY COST 100)

if (${ACTIVATE_RESTART_TESTS})
    add_test(NAME TestSimulationLandauRestartFemUniform_xperiod_vx
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_landau_restart.sh"
            "${PROJECT_SOURCE_DIR}"
            "$<TARGET_FILE:landau_fem_uniform_xperiod_vx>"
            "$<TARGET_FILE:Python3::Interpreter>"
            "restart"
            "${RELATIVE_RESTART_TOLERANCE}"
            "${ABSOLUTE_RESTART_TOLERANCE}")
    set_property(TEST TestSimulationLandauRestartFemUniform_xperiod_vx PROPERTY TIMEOUT 200)
    set_property(TEST TestSimulationLandauRestartFemUniform_xperiod_vx PROPERTY COST 100)
endif()
