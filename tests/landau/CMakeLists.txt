cmake_minimum_required(VERSION 3.15)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(NAME TestLandau
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_landau.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:voicexx>"
        "$<TARGET_FILE:Python3::Interpreter>")
set_property(TEST TestLandau PROPERTY TIMEOUT 200)
set_property(TEST TestLandau PROPERTY LABELS test_on_Release_only)
