# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(NAME TestSimulationLandauFFT_XYVxVy
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_landau4d_small.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:landau4d_fft>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "fft")
set_property(TEST TestSimulationLandauFFT_XYVxVy PROPERTY TIMEOUT 200)
