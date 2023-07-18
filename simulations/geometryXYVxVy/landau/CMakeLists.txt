# SPDX-License-Identifier: MIT

add_executable(landau4d_fft landau4d_fft.cpp)
target_compile_features(landau4d_fft PUBLIC cxx_std_17)
target_link_libraries(landau4d_fft
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
        paraconf::paraconf
        PDI::pdi
        vcx::paraconfpp
        vcx::geometry_xyvxvy
        vcx::initialization_xyvxvy
        vcx::interpolation
        vcx::advection
        vcx::vlasov_xyvxvy
        vcx::poisson_xy
        vcx::time_integration_xyvxvy
        vcx::utils
)

install(TARGETS landau4d_fft)
