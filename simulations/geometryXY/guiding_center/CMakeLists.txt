# SPDX-License-Identifier: MIT

add_executable(guiding_center_XY guiding_center.cpp)
target_link_libraries(guiding_center_XY
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
        paraconf::paraconf
        PDI::pdi

        gslx::advection
        gslx::initialization_Kelvin_Helmholtz
        gslx::interpolation
        gslx::io
        gslx::geometry_XY
        gslx::paraconfpp
        gslx::pde_solvers
        gslx::predcorr_rk2_xy
        gslx::simulation_utils
        gslx::timestepper
        gslx::utils
    
)

install(TARGETS guiding_center_XY)
