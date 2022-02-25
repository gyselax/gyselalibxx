# SPDX-License-Identifier: MIT

add_library(interpolation STATIC
    spline_interpolator_x.cpp
    spline_interpolator_vx.cpp
)

target_compile_features(interpolation
    PUBLIC
        cxx_std_17
)

target_include_directories(interpolation
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(interpolation
    PUBLIC
        DDC::DDC
        sll::splines
        vcx::geometry
)

add_library(vcx::interpolation ALIAS interpolation)
