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
    gslx::geometry_xperiod_vx
)
target_compile_definitions(geometry_xyvxvy INTERFACE GYSELALIBXX_BUILD_REAL_PRECISION=float)
add_library("gslx::geometry_xyvxvy" ALIAS "geometry_xyvxvy")
