# SPDX-License-Identifier: MIT


add_library("geometry_RTheta" INTERFACE)
target_include_directories("geometry_RTheta"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("geometry_RTheta" INTERFACE
    DDC::DDC
    sll::SLL
    gslx::data_types
    gslx::speciesinfo
    gslx::utils
)
add_library("gslx::geometry_RTheta" ALIAS "geometry_RTheta")

