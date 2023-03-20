# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("advection_${GEOMETRY_VARIANT}" STATIC
    bsl_advection_x.cpp
    bsl_advection_vx.cpp
    nulladvectionvx.cpp
)

target_compile_features("advection_${GEOMETRY_VARIANT}"
    PUBLIC
        cxx_std_17
)

target_include_directories("advection_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("advection_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        sll::splines
        vcx::interpolation_${GEOMETRY_VARIANT}
        vcx::speciesinfo
        vcx::geometry_${GEOMETRY_VARIANT}
)

add_library("vcx::advection_${GEOMETRY_VARIANT}" ALIAS "advection_${GEOMETRY_VARIANT}")

endforeach()
