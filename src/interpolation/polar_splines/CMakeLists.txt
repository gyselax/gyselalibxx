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
        gyselalibxx::coord_transformations
        gyselalibxx::math_tools
        gyselalibxx::utils
)

add_library("gyselalibxx::polar_splines" ALIAS "polar_splines")

install(TARGETS polar_splines EXPORT gyselalibxxTargets FILE_SET HEADERS)
