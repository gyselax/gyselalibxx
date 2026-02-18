# SPDX-License-Identifier: MIT

set(BASIC_GEOMETRY_XVx_VARIANTS_LIST xperiod_vx xnonperiod_vx)

foreach(GEOMETRY_VARIANT IN LISTS BASIC_GEOMETRY_XVx_VARIANTS_LIST)

add_executable(vlasovpoisson_xvx_fem_uniform_${GEOMETRY_VARIANT} vlasovpoisson_xvx_fem_uniform.cpp)
target_link_libraries(vlasovpoisson_xvx_fem_uniform_${GEOMETRY_VARIANT}
    PUBLIC
        DDC::core
        DDC::pdi
        paraconf::paraconf
        PDI::pdi
        gslx::initialisation_${GEOMETRY_VARIANT}
        gslx::interpolation
        gslx::paraconfpp
        gslx::poisson_${GEOMETRY_VARIANT}
        gslx::speciesinfo
        gslx::time_integration_${GEOMETRY_VARIANT}
        gslx::boltzmann_${GEOMETRY_VARIANT}
        gslx::advection
        gslx::io
        gslx::pde_solvers
        gslx::utils

)

install(TARGETS vlasovpoisson_xvx_fem_uniform_${GEOMETRY_VARIANT})

endforeach()

add_executable(vlasovpoisson_xvx_fft vlasovpoisson_xvx_fft.cpp)
target_link_libraries(vlasovpoisson_xvx_fft
    PUBLIC
        DDC::core
        DDC::pdi
        paraconf::paraconf
        PDI::pdi
        gslx::initialisation_xperiod_vx
        gslx::interpolation
        gslx::paraconfpp
        gslx::poisson_xperiod_vx
        gslx::speciesinfo
        gslx::time_integration_xperiod_vx
        gslx::boltzmann_xperiod_vx
        gslx::advection
        gslx::io
        gslx::pde_solvers
        gslx::utils

)

install(TARGETS vlasovpoisson_xvx_fft)
