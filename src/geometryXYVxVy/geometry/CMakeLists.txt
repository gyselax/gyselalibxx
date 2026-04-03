# SPDX-License-Identifier: MIT

foreach(PRECISION_VARIANT IN LISTS PRECISION_VARIANTS_LIST)

add_library("geometry_xyvxvy_${PRECISION_VARIANT}" INTERFACE)

target_compile_features("geometry_xyvxvy_${PRECISION_VARIANT}" INTERFACE cxx_std_20)

target_sources("geometry_xyvxvy_${PRECISION_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                geometry_xyvxvy.hpp
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

install(TARGETS geometry_xyvxvy_${PRECISION_VARIANT} EXPORT GyselalibxxTargets FILE_SET HEADERS)

endforeach()
