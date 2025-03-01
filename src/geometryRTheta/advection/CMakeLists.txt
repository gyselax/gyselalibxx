# SPDX-License-Identifier: MIT

add_library(advection_polar INTERFACE)
target_include_directories(advection_polar
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(advection_polar 
    INTERFACE
        DDC::core
        Eigen3::Eigen
        gslx::advection
        gslx::data_types
        gslx::geometry_RTheta
        gslx::interpolation_polar
        gslx::mapping
        gslx::math_tools
        gslx::utils
)
add_library(gslx::advection_polar ALIAS advection_polar)
