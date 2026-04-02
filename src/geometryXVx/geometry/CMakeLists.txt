# SPDX-License-Identifier: MIT


foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("geometry_${GEOMETRY_VARIANT}" INTERFACE)

target_sources("geometry_${GEOMETRY_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                geometry_xvx.hpp
)

target_link_libraries("geometry_${GEOMETRY_VARIANT}" INTERFACE
    DDC::core
    gslx::speciesinfo
    gslx::moments
    gslx::utils
)
add_library("gslx::geometry_${GEOMETRY_VARIANT}" ALIAS "geometry_${GEOMETRY_VARIANT}")

install(TARGETS geometry_${GEOMETRY_VARIANT} EXPORT GyselalibxxTargets FILE_SET HEADERS)

endforeach()

target_compile_definitions(geometry_xperiod_vx INTERFACE PERIODIC_RDIMX)
target_compile_definitions(geometry_xperiod_vx_user_defined_grid INTERFACE PERIODIC_RDIMX)
target_compile_definitions(geometry_xperiod_vx_user_defined_grid INTERFACE INPUT_MESH)
target_compile_definitions(geometry_xnonperiod_vx_user_defined_grid INTERFACE INPUT_MESH)
