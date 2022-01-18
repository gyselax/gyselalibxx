add_executable(landau_fft
    landau_fft.cpp
)
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
        FFTW::Double
        PDI::pdi
        sll::splines
)

add_executable(landau_fem_periodic
    landau_fem_periodic.cpp
)
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
        FFTW::Double
        PDI::pdi
        sll::splines
)
