# SPDX-License-Identifier: MIT


add_library(time_integration_rtheta INTERFACE)
target_include_directories(time_integration_rtheta
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(time_integration_rtheta 
    INTERFACE
        DDC::core
        Eigen3::Eigen
        gslx::interpolation_RTheta
        gslx::geometry_RTheta
        gslx::poisson_RTheta
        gslx::advection_RTheta
        gslx::utils
        gslx::timestepper
        gslx::advection_field_RTheta
)
add_library(gslx::time_integration_rtheta ALIAS time_integration_rtheta)

            
