# SPDX-License-Identifier: MIT

set(BASIC_GEOMETRY_XVx_VARIANTS_LIST xperiod_vx xnonperiod_vx)

foreach(GEOMETRY_VARIANT IN LISTS BASIC_GEOMETRY_XVx_VARIANTS_LIST)

add_executable(bumpontail_fem_uniform_${GEOMETRY_VARIANT} bumpontail_fem_uniform.cpp)
target_link_libraries(bumpontail_fem_uniform_${GEOMETRY_VARIANT}
    PUBLIC
        DDC::core
        DDC::pdi
        paraconf::paraconf
        PDI::pdi
        gslx::initialization_${GEOMETRY_VARIANT}
        gslx::interpolation
        gslx::paraconfpp
        gslx::poisson_${GEOMETRY_VARIANT}
        gslx::speciesinfo
        gslx::time_integration_${GEOMETRY_VARIANT}
        gslx::boltzmann_${GEOMETRY_VARIANT}
        gslx::advection
        gslx::io
        gslx::utils

)

install(TARGETS bumpontail_fem_uniform_${GEOMETRY_VARIANT})

endforeach()

add_executable(bumpontail_fft bumpontail_fft.cpp)
target_link_libraries(bumpontail_fft
    PUBLIC
        DDC::core
        DDC::pdi
        paraconf::paraconf
        PDI::pdi
        gslx::initialization_xperiod_vx
        gslx::interpolation
        gslx::paraconfpp
        gslx::poisson_xperiod_vx
        gslx::speciesinfo
        gslx::time_integration_xperiod_vx
        gslx::boltzmann_xperiod_vx
        gslx::advection
        gslx::io
        gslx::utils

)

install(TARGETS bumpontail_fft)
