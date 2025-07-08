# SPDX-License-Identifier: MIT

add_library(testing_tools STATIC
    test_math_utils.cpp
)

target_include_directories(testing_tools
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)
