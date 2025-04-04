# SPDX-License-Identifier: MIT


add_library(time_integration_RTheta INTERFACE)
target_include_directories(time_integration_RTheta
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(time_integration_RTheta 
    INTERFACE
        DDC::core
        Eigen3::Eigen
        gslx::advection
        gslx::advection_field_RTheta
        gslx::geometry_RTheta
        gslx::poisson_RTheta
        gslx::utils
        gslx::timestepper
)
add_library(gslx::time_integration_RTheta ALIAS time_integration_RTheta)

            
