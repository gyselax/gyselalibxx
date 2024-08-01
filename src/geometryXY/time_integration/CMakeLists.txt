# SPDX-License-Identifier: MIT

add_library("predcorr_rk2_xy" INTERFACE)
target_include_directories("predcorr_rk2_xy"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("predcorr_rk2_xy" INTERFACE
    DDC::DDC
    
    gslx::advection
    gslx::interpolation
    gslx::geometry_XY
    gslx::pde_solvers
    gslx::timestepper
    gslx::utils



)
add_library("gslx::predcorr_rk2_xy" ALIAS "predcorr_rk2_xy")
