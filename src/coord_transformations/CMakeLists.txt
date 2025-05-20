# SPDX-License-Identifier: MIT

add_library("coord_transformations" INTERFACE)

target_include_directories("coord_transformations"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("coord_transformations"
    INTERFACE
        DDC::core
        gslx::data_types
        gslx::math_tools
        gslx::utils
)

add_library("gslx::coord_transformations" ALIAS "coord_transformations")

