# SPDX-License-Identifier: MIT

add_executable(guiding_centre_XY guiding_centre.cpp)
target_link_libraries(guiding_centre_XY
    PUBLIC
        DDC::core
        DDC::pdi
        paraconf::paraconf
        PDI::PDI_C

        gyselalibxx::advection
        gyselalibxx::initialisation_Kelvin_Helmholtz
        gyselalibxx::interpolation
        gyselalibxx::io
        gyselalibxx::geometry_XY
        gyselalibxx::paraconfpp
        gyselalibxx::pde_solvers
        gyselalibxx::predcorr_rk2_xy
        gyselalibxx::simulation_utils
        gyselalibxx::timestepper
        gyselalibxx::utils
    
)

install(TARGETS guiding_centre_XY)
