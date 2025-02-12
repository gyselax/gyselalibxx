# SPDX-License-Identifier: MIT

add_library(speciesinfo INTERFACE)
target_include_directories(speciesinfo
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_link_libraries(speciesinfo INTERFACE
    DDC::core
    gslx::utils
)
add_library(gslx::speciesinfo ALIAS speciesinfo)


add_library(moments INTERFACE)
target_include_directories(moments
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_link_libraries(moments INTERFACE
    DDC::core)
add_library(gslx::moments ALIAS moments)
