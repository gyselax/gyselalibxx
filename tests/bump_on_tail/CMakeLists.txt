# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

if(VOICEXX_ENABLE_PERIODIC_RDIMX)
    add_test(NAME TestBumpontailFFT
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_bumpontail.sh"
            "${PROJECT_SOURCE_DIR}"
            "$<TARGET_FILE:bumpontail_fft>"
            "$<TARGET_FILE:Python3::Interpreter>"
            "fft")
    set_property(TEST TestBumpontailFFT PROPERTY TIMEOUT 200)
    set_property(TEST TestBumpontailFFT PROPERTY LABELS test_on_Release_only)
endif()

add_test(NAME TestBumpontailFemUniform
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_bumpontail.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:bumpontail_fem_uniform>"
        "$<TARGET_FILE:Python3::Interpreter>"
        "fem")
set_property(TEST TestBumpontailFemUniform PROPERTY TIMEOUT 200)
set_property(TEST TestBumpontailFemUniform PROPERTY LABELS test_on_Release_only)
