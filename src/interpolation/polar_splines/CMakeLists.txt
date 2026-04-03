# SPDX-License-Identifier: MIT

add_library("polar_splines" INTERFACE)

target_compile_features("polar_splines" INTERFACE cxx_std_20)

target_sources("polar_splines"
    INTERFACE
        FILE_SET HEADERS
            FILES
                bernstein.hpp
                polar_bsplines.hpp
                polar_spline_evaluator.hpp
)

target_link_libraries("polar_splines"
    INTERFACE
        DDC::core
        DDC::splines
        gslx::coord_transformations
        gslx::math_tools
        gslx::utils
)

add_library("gslx::polar_splines" ALIAS "polar_splines")

install(TARGETS polar_splines EXPORT GyselalibxxTargets FILE_SET HEADERS)
