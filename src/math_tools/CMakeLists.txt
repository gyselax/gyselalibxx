# SPDX-License-Identifier: MIT

add_library("math_tools" INTERFACE)

target_compile_features("math_tools" INTERFACE cxx_std_20)

target_sources("math_tools"
    INTERFACE
        FILE_SET HEADERS
            FILES
                central_fdm_partial_derivatives.hpp
                central_fdm_partial_derivatives_with_boundary_values.hpp
                constant_partial_derivatives.hpp
                curl.hpp
                gradient_creator.hpp
                gradient.hpp
                ipartial_derivative.hpp
                lie_poisson_bracket.hpp
                l_norm_tools.hpp
                math_tools.hpp
                spline_1d_partial_derivative.hpp
                spline_2d_partial_derivative.hpp
                spline_builder_2d_cache.hpp
)

target_link_libraries("math_tools"
    INTERFACE
        DDC::core
        gslx::coord_transformations
        gslx::data_types
        gslx::interpolation
        gslx::quadrature
)

add_library("gslx::math_tools" ALIAS "math_tools")

install(TARGETS math_tools EXPORT GyselalibxxTargets FILE_SET HEADERS)
