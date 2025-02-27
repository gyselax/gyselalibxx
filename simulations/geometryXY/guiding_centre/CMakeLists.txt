# SPDX-License-Identifier: MIT

add_executable(guiding_centre_XY guiding_centre.cpp)
target_link_libraries(guiding_centre_XY
    PUBLIC
        DDC::core
        DDC::pdi
        paraconf::paraconf
        PDI::pdi

        gslx::advection
        gslx::initialisation_Kelvin_Helmholtz
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

install(TARGETS guiding_centre_XY)
