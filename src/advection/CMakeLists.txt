# SPDX-License-Identifier: MIT

add_library("advection" INTERFACE)

target_compile_features("advection" INTERFACE cxx_std_20)

target_sources("advection"
    INTERFACE
        FILE_SET HEADERS
            FILES
                bsl_advection_1d.hpp
                bsl_advection_polar.hpp
                bsl_advection_vx.hpp
                bsl_advection_x.hpp
                iadvectionvx.hpp
                iadvectionx.hpp
                ipolar_foot_finder.hpp
                nulladvectionvx.hpp
                spline_polar_foot_finder.hpp
)

target_link_libraries("advection"
    INTERFACE
        DDC::core
        gslx::interpolation
        gslx::math_tools
        gslx::speciesinfo
        gslx::timestepper
        gslx::utils
)

add_library("gslx::advection" ALIAS "advection")

install(TARGETS advection EXPORT GyselalibxxTargets FILE_SET HEADERS)
