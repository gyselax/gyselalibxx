# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("interpolation_${GEOMETRY_VARIANT}" STATIC
    spline_interpolator_x.cpp
    spline_interpolator_vx.cpp
)

target_compile_features("interpolation_${GEOMETRY_VARIANT}"
    PUBLIC
        cxx_std_17
)

target_include_directories("interpolation_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("interpolation_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        sll::splines
        vcx::geometry_${GEOMETRY_VARIANT}
)

add_library("vcx::interpolation_${GEOMETRY_VARIANT}" ALIAS "interpolation_${GEOMETRY_VARIANT}")

endforeach()
