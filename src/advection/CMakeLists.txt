# SPDX-License-Identifier: MIT

add_library("advection" INTERFACE)

target_compile_features("advection" INTERFACE cxx_std_20)

target_sources("advection"
    INTERFACE
        FILE_SET HEADERS
            FILES
                polar_foot_finders/elementwise_choice.hpp
                polar_foot_finders/logical_advection_logical_foot_finder.hpp
                polar_foot_finders/logical_advection_pseudo_physical_foot_finder.hpp
                polar_foot_finders/physical_advection_physical_foot_finder.hpp
                polar_foot_finders/physical_advection_pseudo_physical_foot_finder.hpp
                bsl_advection_1d.hpp
                bsl_advection_polar.hpp
                bsl_advection_vx.hpp
                bsl_advection_x.hpp
                iadvectionvx.hpp
                iadvectionx.hpp
                nulladvectionvx.hpp
                polar_foot_finder.hpp
                spline_polar_foot_finder.hpp
)

target_link_libraries("advection"
    INTERFACE
        DDC::core
        gyselalibxx::interpolation
        gyselalibxx::math_tools
        gyselalibxx::speciesinfo
        gyselalibxx::timestepper
        gyselalibxx::utils
)

add_library("gyselalibxx::advection" ALIAS "advection")

install(TARGETS advection EXPORT GyselalibxxTargets FILE_SET HEADERS)
