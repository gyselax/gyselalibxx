# SPDX-License-Identifier: MIT

add_library(advection_rtheta INTERFACE)
target_include_directories(advection_rtheta
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(advection_rtheta 
    INTERFACE
        DDC::core
        Eigen3::Eigen
        gslx::advection
        gslx::data_types
        gslx::geometry_RTheta
        gslx::interpolation_rtheta
        gslx::mapping
        gslx::math_tools
        gslx::utils
)
add_library(gslx::advection_rtheta ALIAS advection_rtheta)
