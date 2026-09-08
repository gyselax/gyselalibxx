# SPDX-License-Identifier: MIT

set(BASIC_GEOMETRY_XVx_VARIANTS_LIST xperiod_vx xnonperiod_vx)

foreach(GEOMETRY_VARIANT IN LISTS BASIC_GEOMETRY_XVx_VARIANTS_LIST)

add_executable(vlasovpoisson_xvx_fem_uniform_${GEOMETRY_VARIANT} vlasovpoisson_xvx_fem_uniform.cpp)
target_link_libraries(vlasovpoisson_xvx_fem_uniform_${GEOMETRY_VARIANT}
    PUBLIC
        DDC::core
        DDC::pdi
        paraconf::paraconf
        PDI::PDI_C
        gyselalibxx::initialisation_${GEOMETRY_VARIANT}
        gyselalibxx::interpolation
        gyselalibxx::paraconfpp
        gyselalibxx::poisson_${GEOMETRY_VARIANT}
        gyselalibxx::speciesinfo
        gyselalibxx::time_integration_${GEOMETRY_VARIANT}
        gyselalibxx::boltzmann_${GEOMETRY_VARIANT}
        gyselalibxx::advection
        gyselalibxx::io
        gyselalibxx::pde_solvers
        gyselalibxx::utils

)

install(TARGETS vlasovpoisson_xvx_fem_uniform_${GEOMETRY_VARIANT})

endforeach()

add_executable(vlasovpoisson_xvx_fft vlasovpoisson_xvx_fft.cpp)
target_link_libraries(vlasovpoisson_xvx_fft
    PUBLIC
        DDC::core
        DDC::pdi
        paraconf::paraconf
        PDI::PDI_C
        gyselalibxx::initialisation_xperiod_vx
        gyselalibxx::interpolation
        gyselalibxx::paraconfpp
        gyselalibxx::poisson_xperiod_vx
        gyselalibxx::speciesinfo
        gyselalibxx::time_integration_xperiod_vx
        gyselalibxx::boltzmann_xperiod_vx
        gyselalibxx::advection
        gyselalibxx::io
        gyselalibxx::pde_solvers
        gyselalibxx::utils

)

install(TARGETS vlasovpoisson_xvx_fft)
