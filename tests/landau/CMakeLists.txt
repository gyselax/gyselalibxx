# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(NAME TestLandauFFT
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_landau.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:landau_fft>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "fft")
set_property(TEST TestLandauFFT PROPERTY TIMEOUT 200)

add_test(NAME TestLandauFemUniform_periodic_rdimx
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_landau.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:landau_fem_uniform_periodic_rdimx>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "fem")
set_property(TEST TestLandauFemUniform_periodic_rdimx PROPERTY TIMEOUT 200)
