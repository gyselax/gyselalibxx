# SPDX-License-Identifier: MIT

add_executable(bumpontail_fem_uniform bumpontail_fem_uniform.cpp)
target_compile_features(bumpontail_fem_uniform PUBLIC cxx_std_17)
target_link_libraries(bumpontail_fem_uniform
  PUBLIC
      DDC::DDC
      DDC::PDI_Wrapper
      vcx::initialization
      vcx::interpolation
      vcx::poisson
      vcx::time
      vcx::speciesinfo
      vcx::vlasov
      vcx::paraconfpp
      PDI::pdi
      paraconf::paraconf
      sll::splines
)

install(TARGETS bumpontail_fem_uniform)

if (VOICEXX_ENABLE_PERIODIC_RDIMX)
    add_executable(bumpontail_fft bumpontail_fft.cpp)
    target_compile_features(bumpontail_fft PUBLIC cxx_std_17)
    target_link_libraries(bumpontail_fft
        PUBLIC
            DDC::DDC
            DDC::PDI_Wrapper
            vcx::initialization
            vcx::fft
            vcx::interpolation
            vcx::poisson
            vcx::time
            vcx::speciesinfo
            vcx::vlasov
            vcx::paraconfpp
            PDI::pdi
            paraconf::paraconf
            sll::splines
    )

    install(TARGETS bumpontail_fft)
endif()
