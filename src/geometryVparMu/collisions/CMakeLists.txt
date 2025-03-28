# SPDX-License-Identifier: MIT

add_library("geometry_collisions_vparmu" INTERFACE)

target_include_directories("geometry_collisions_vparmu"
    INTERFACE
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("geometry_collisions_vparmu"
    INTERFACE
        gslx::collisions
)

add_library("gslx::geometry_collisions_vparmu" ALIAS "geometry_collisions_vparmu")
