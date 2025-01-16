# SPDX-License-Identifier: MIT

add_library(advection_rp INTERFACE)
target_include_directories(advection_rp
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(advection_rp 
    INTERFACE
        DDC::core
        Eigen3::Eigen
        gslx::data_types
        gslx::geometry_RTheta
        gslx::interpolation_2D_rp
        gslx::mapping
        gslx::math_tools
        gslx::utils
)
add_library(gslx::advection_rp ALIAS advection_rp)
