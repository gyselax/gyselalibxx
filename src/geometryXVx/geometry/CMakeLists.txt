# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("geometry_${GEOMETRY_VARIANT}" INTERFACE)
target_include_directories("geometry_${GEOMETRY_VARIANT}"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("geometry_${GEOMETRY_VARIANT}" INTERFACE
    DDC::DDC
    sll::splines
    gslx::speciesinfo
)
add_library("gslx::geometry_${GEOMETRY_VARIANT}" ALIAS "geometry_${GEOMETRY_VARIANT}")

endforeach()

target_compile_definitions(geometry_xperiod_vx INTERFACE PERIODIC_RDIMX)
