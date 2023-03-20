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
        sll::splines
        vcx::initialization_${GEOMETRY_VARIANT}
        vcx::interpolation_${GEOMETRY_VARIANT}
        vcx::paraconfpp
        vcx::poisson_${GEOMETRY_VARIANT}
        vcx::speciesinfo
        vcx::time_integration_${GEOMETRY_VARIANT}
        vcx::boltzmann_${GEOMETRY_VARIANT}
        vcx::advection_${GEOMETRY_VARIANT}
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
        sll::splines
        vcx::fft
        vcx::initialization_xperiod_vx
        vcx::interpolation_xperiod_vx
        vcx::paraconfpp
        vcx::poisson_xperiod_vx
        vcx::speciesinfo
        vcx::time_integration_xperiod_vx
        vcx::boltzmann_xperiod_vx
        vcx::advection_xperiod_vx
)

install(TARGETS landau_fft)
