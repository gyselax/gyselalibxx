# SPDX-License-Identifier: MIT


add_library("geometry_xvxvyvz" INTERFACE)
target_include_directories("geometry_xvxvyvz"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("geometry_xvxvyvz" INTERFACE
    DDC::core
    gslx::mpi_parallelisation
    gslx::speciesinfo
    gslx::utils
)
add_library("gslx::geometry_xvxvyvz" ALIAS "geometry_xvxvyvz")
