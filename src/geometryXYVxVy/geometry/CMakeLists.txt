# SPDX-License-Identifier: MIT

foreach(PRECISION_VARIANT IN LISTS PRECISION_VARIANTS_LIST)

add_library("geometry_xyvxvy_${PRECISION_VARIANT}" INTERFACE)
target_include_directories("geometry_xyvxvy_${PRECISION_VARIANT}"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("geometry_xyvxvy_${PRECISION_VARIANT}" INTERFACE
    DDC::core
    gslx::data_types
    gslx::mpi_parallelisation
    gslx::speciesinfo
    gslx::utils
    gslx::geometry_xperiod_vx
)
target_compile_definitions(geometry_xyvxvy_${PRECISION_VARIANT} INTERFACE GYSELALIBXX_BUILD_REAL_PRECISION=${PRECISION_VARIANT})
add_library("gslx::geometry_xyvxvy_${PRECISION_VARIANT}" ALIAS "geometry_xyvxvy_${PRECISION_VARIANT}")

endforeach()
