# SPDX-License-Identifier: MIT

add_library(speciesinfo STATIC
    species_init.cpp
)
target_include_directories(speciesinfo
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_link_libraries(speciesinfo 
    DDC::DDC
    gslx::io
    gslx::paraconfpp
    gslx::utils
)
add_library(gslx::speciesinfo ALIAS speciesinfo)


add_library(moments INTERFACE)
target_include_directories(moments
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_link_libraries(moments INTERFACE
    DDC::DDC)
add_library(gslx::moments ALIAS moments)
