# SPDX-License-Identifier: MIT


add_library("geometry_xyvxvy" INTERFACE)
target_include_directories("geometry_xyvxvy"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("geometry_xyvxvy" INTERFACE
    DDC::core
    gslx::data_types
    gslx::mpi_parallelisation
    gslx::speciesinfo
    gslx::utils
)
add_library("gslx::geometry_xyvxvy" ALIAS "geometry_xyvxvy")
