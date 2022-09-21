# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_VARIANTS_LIST)

add_executable(bumpontail_fem_uniform_${GEOMETRY_VARIANT} bumpontail_fem_uniform.cpp)
target_compile_features(bumpontail_fem_uniform_${GEOMETRY_VARIANT} PUBLIC cxx_std_17)
target_link_libraries(bumpontail_fem_uniform_${GEOMETRY_VARIANT}
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
        vcx::time_${GEOMETRY_VARIANT}
        vcx::boltzmann_${GEOMETRY_VARIANT}
)

install(TARGETS bumpontail_fem_uniform_${GEOMETRY_VARIANT})

endforeach()

add_executable(bumpontail_fft bumpontail_fft.cpp)
target_compile_features(bumpontail_fft PUBLIC cxx_std_17)
target_link_libraries(bumpontail_fft
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
        paraconf::paraconf
        PDI::pdi
        sll::splines
        vcx::fft
        vcx::initialization_periodic_rdimx
        vcx::interpolation_periodic_rdimx
        vcx::paraconfpp
        vcx::poisson_periodic_rdimx
        vcx::speciesinfo
        vcx::time_periodic_rdimx
        vcx::boltzmann_periodic_rdimx
)

install(TARGETS bumpontail_fft)
