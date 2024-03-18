# SPDX-License-Identifier: MIT

add_executable(landau4d_fft landau4d_fft.cpp)
target_compile_features(landau4d_fft PUBLIC cxx_std_17)
target_link_libraries(landau4d_fft
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
        paraconf::paraconf
        PDI::pdi
        gslx::paraconfpp
        gslx::geometry_xyvxvy
        gslx::initialization_xyvxvy
        gslx::interpolation
        gslx::advection
        gslx::vlasov_xyvxvy
        gslx::poisson_xy
        gslx::quadrature
        gslx::time_integration_xyvxvy
        gslx::utils
)

install(TARGETS landau4d_fft)
