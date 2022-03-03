# SPDX-License-Identifier: MIT

add_executable(landau_fem_uniform landau_fem_uniform.cpp)
target_compile_features(landau_fem_uniform PUBLIC cxx_std_17)
target_link_libraries(landau_fem_uniform
  PUBLIC
  DDC::DDC
  DDC::PDI
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

install(TARGETS landau_fem_uniform)

if (VOICEXX_ENABLE_PERIODIC_RDIMX)
    add_executable(landau_fft landau_fft.cpp)
    target_compile_features(landau_fft PUBLIC cxx_std_17)
    target_link_libraries(landau_fft
        PUBLIC
            DDC::DDC
            DDC::PDI
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

    install(TARGETS landau_fft)
endif()
