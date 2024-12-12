# SPDX-License-Identifier: MIT

add_library("geometry_XY" INTERFACE)
target_include_directories("geometry_XY"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("geometry_XY" INTERFACE
    DDC::core
    DDC::splines
    gslx::utils
)
add_library("gslx::geometry_XY" ALIAS "geometry_XY")
