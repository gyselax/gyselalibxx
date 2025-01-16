# SPDX-License-Identifier: MIT

add_library("geometry_vparmu" INTERFACE)
target_include_directories("geometry_vparmu"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("geometry_vparmu" INTERFACE
    DDC::core
    gslx::speciesinfo
    gslx::utils
)
add_library("gslx::geometry_vparmu" ALIAS "geometry_vparmu")

