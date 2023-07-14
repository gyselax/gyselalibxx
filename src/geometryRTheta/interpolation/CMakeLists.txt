# SPDX-License-Identifier: MIT


add_library("interpolation_2D_rp" STATIC
    spline_interpolator_2d_rp.cpp
)

target_compile_features("interpolation_2D_rp"
    PUBLIC
        cxx_std_17
)

target_include_directories("interpolation_2D_rp"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("interpolation_2D_rp"
    PUBLIC
        DDC::DDC
        sll::splines
        vcx::geometry_RTheta
)

add_library("vcx::interpolation_2D_rp" ALIAS "interpolation_2D_rp")
