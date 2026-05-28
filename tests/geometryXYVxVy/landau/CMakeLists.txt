# SPDX-License-Identifier: MIT


find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(NAME TestSimulationLandauFFT_XYVxVy_spline
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_landau4d_small.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:landau4d_fft_spline>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "fft"
        "${CMAKE_CURRENT_SOURCE_DIR}/landau_small.yaml")
set_property(TEST TestSimulationLandauFFT_XYVxVy_spline PROPERTY TIMEOUT 200)
set_property(TEST TestSimulationLandauFFT_XYVxVy_spline PROPERTY COST 100)

add_test(NAME TestSimulationLandauFFT_XYVxVy_lagrange
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_landau4d_small.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:landau4d_fft_lagrange>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "fft"
        "${CMAKE_CURRENT_SOURCE_DIR}/landau_small_lagrange.yaml")
set_property(TEST TestSimulationLandauFFT_XYVxVy_lagrange PROPERTY TIMEOUT 200)
set_property(TEST TestSimulationLandauFFT_XYVxVy_lagrange PROPERTY COST 100)
