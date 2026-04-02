# SPDX-License-Identifier: MIT

add_library("interpolation" INTERFACE)

target_sources("interpolation"
    INTERFACE
        FILE_SET HEADERS
            FILES
                constant_identity_interpolation_extrapolation_rule.hpp
                extrapolation_rule_choice.hpp
                identity_interpolation_builder.hpp
                i_interpolation_builder.hpp
                i_interpolation_evaluator.hpp
                i_interpolation.hpp
                lagrange_basis_non_uniform.hpp
                lagrange_basis_uniform.hpp
                lagrange_evaluator.hpp
                lagrange_interpolation.hpp
                nd_lagrange_evaluator.hpp
                spline_interpolation.hpp
)

target_link_libraries("interpolation"
    INTERFACE
        DDC::core
        gslx::math_tools
        gslx::utils
)

add_library("gslx::interpolation" ALIAS "interpolation")

install(TARGETS interpolation EXPORT GyselalibxxTargets FILE_SET HEADERS)

add_subdirectory(polar_splines)
