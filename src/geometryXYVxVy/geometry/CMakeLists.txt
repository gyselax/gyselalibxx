# SPDX-License-Identifier: MIT


add_library("geometry_xyvxvy" INTERFACE)

target_sources("geometry_xyvxvy"
    INTERFACE
        FILE_SET HEADERS
            FILES
                geometry_xyvxvy.hpp
)

target_link_libraries("geometry_xyvxvy" INTERFACE
    DDC::core
    gslx::data_types
    gslx::mpi_parallelisation
    gslx::speciesinfo
    gslx::utils
    gslx::geometry_xperiod_vx
)
add_library("gslx::geometry_xyvxvy" ALIAS "geometry_xyvxvy")

install(TARGETS geometry_xyvxvy EXPORT GyselalibxxTargets FILE_SET HEADERS)
