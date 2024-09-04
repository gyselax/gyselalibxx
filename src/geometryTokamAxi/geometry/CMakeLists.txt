# SPDX-License-Identifier: MIT

add_library("geometry_tokamaxi" INTERFACE)
target_include_directories("geometry_tokamaxi"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("geometry_tokamaxi" INTERFACE
    DDC::DDC
    gslx::speciesinfo
)
add_library("gslx::geometry_tokamaxi" ALIAS "geometry_tokamaxi")

