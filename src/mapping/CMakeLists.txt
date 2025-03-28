# SPDX-License-Identifier: MIT

add_library("mapping" INTERFACE)

target_include_directories("mapping"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("mapping"
    INTERFACE
        DDC::core
        gslx::data_types
        gslx::math_tools
        gslx::utils
)

add_library("gslx::mapping" ALIAS "mapping")

