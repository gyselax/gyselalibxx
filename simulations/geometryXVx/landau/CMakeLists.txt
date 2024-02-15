# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_executable(landau_fem_uniform_${GEOMETRY_VARIANT} landau_fem_uniform.cpp)
target_compile_features(landau_fem_uniform_${GEOMETRY_VARIANT} PUBLIC cxx_std_17)
target_link_libraries(landau_fem_uniform_${GEOMETRY_VARIANT}
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
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
)

install(TARGETS landau_fem_uniform_${GEOMETRY_VARIANT})

endforeach()

add_executable(landau_fft landau_fft.cpp)
target_compile_features(landau_fft PUBLIC cxx_std_17)
target_link_libraries(landau_fft
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
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
)

install(TARGETS landau_fft)
