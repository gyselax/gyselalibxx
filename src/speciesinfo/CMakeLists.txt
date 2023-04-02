# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

add_library(speciesinfo INTERFACE)
target_include_directories(speciesinfo
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(speciesinfo INTERFACE
    DDC::DDC
    sll::splines
)
add_library(vcx::speciesinfo ALIAS speciesinfo)
