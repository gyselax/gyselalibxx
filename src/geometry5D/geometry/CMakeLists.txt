# SPDX-License-Identifier: MIT

add_library("geometry5D" INTERFACE)
target_include_directories("geometry5D"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("geometry5D" INTERFACE
    DDC::DDC
    sll::SLL
    gslx::speciesinfo
    gslx::utils
)
add_library("gslx::geometry5D" ALIAS "geometry5D")

