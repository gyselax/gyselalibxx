# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(NAME TestBumpontailFFT
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_bumpontail.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:bumpontail_fft>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "fft")
set_property(TEST TestBumpontailFFT PROPERTY TIMEOUT 200)

add_test(NAME TestBumpontailFemUniform_periodic_rdimx
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_bumpontail.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:bumpontail_fem_uniform_periodic_rdimx>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "fem")
set_property(TEST TestBumpontailFemUniform_periodic_rdimx PROPERTY TIMEOUT 200)
