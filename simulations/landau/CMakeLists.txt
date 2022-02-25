# SPDX-License-Identifier: MIT

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

    add_executable(landau_fem_periodic landau_fem_periodic.cpp)
    target_compile_features(landau_fem_periodic PUBLIC cxx_std_17)
    target_link_libraries(landau_fem_periodic
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

    install(TARGETS landau_fem_periodic)
else()
    add_executable(landau_fem_non_periodic landau_fem_non_periodic.cpp)
    target_compile_features(landau_fem_non_periodic PUBLIC cxx_std_17)
    target_link_libraries(landau_fem_non_periodic
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
            sll::splines
    )
endif()
