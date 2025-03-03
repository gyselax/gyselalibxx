# SPDX-License-Identifier: MIT

add_library(advection_RTheta INTERFACE)
target_include_directories(advection_RTheta
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(advection_RTheta 
    INTERFACE
        DDC::core
        Eigen3::Eigen
        gslx::advection
        gslx::data_types
        gslx::geometry_RTheta
        gslx::interpolation_RTheta
        gslx::mapping
        gslx::math_tools
        gslx::utils
)
add_library(gslx::advection_RTheta ALIAS advection_RTheta)
