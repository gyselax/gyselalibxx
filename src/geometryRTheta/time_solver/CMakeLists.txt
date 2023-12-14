# SPDX-License-Identifier: MIT


add_library(time_integration_rp INTERFACE)
target_include_directories(time_integration_rp
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(time_integration_rp 
    INTERFACE
        DDC::DDC
        sll::splines
        Eigen3::Eigen
        gslx::interpolation_2D_rp
        gslx::geometry_RTheta
        gslx::poisson_RTheta
        gslx::advection_rp
        gslx::utils
        gslx::timestepper
)
add_library(gslx::time_integration_rp ALIAS time_integration_rp)

            
