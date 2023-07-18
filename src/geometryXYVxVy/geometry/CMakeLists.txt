# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

add_library("geometry_xyvxvy" INTERFACE)
target_include_directories("geometry_xyvxvy"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("geometry_xyvxvy" INTERFACE
    DDC::DDC
    sll::splines
    vcx::speciesinfo
)
add_library("vcx::geometry_xyvxvy" ALIAS "geometry_xyvxvy")
