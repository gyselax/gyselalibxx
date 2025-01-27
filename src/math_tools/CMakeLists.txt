# SPDX-License-Identifier: MIT

add_library("math_tools" INTERFACE)

target_include_directories("math_tools"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("math_tools"
    INTERFACE
        DDC::core
        gslx::data_types
)

add_library("gslx::math_tools" ALIAS "math_tools")

