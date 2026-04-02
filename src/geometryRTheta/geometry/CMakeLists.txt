# SPDX-License-Identifier: MIT


add_library("geometry_RTheta" INTERFACE)

target_sources("geometry_RTheta"
    INTERFACE
        FILE_SET HEADERS
            FILES
                geometry_r_theta.hpp
                spline_definitions_r_theta.hpp
)

target_link_libraries("geometry_RTheta" INTERFACE
    DDC::core
    gslx::data_types
    gslx::polar_splines
    gslx::utils
)
add_library("gslx::geometry_RTheta" ALIAS "geometry_RTheta")

install(TARGETS geometry_RTheta EXPORT GyselalibxxTargets FILE_SET HEADERS)
