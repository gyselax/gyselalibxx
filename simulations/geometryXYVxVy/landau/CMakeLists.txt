# SPDX-License-Identifier: MIT

add_executable(landau4d_fft landau4d_fft.cpp)
target_link_libraries(landau4d_fft
    PUBLIC
        DDC::core
        DDC::pdi
        paraconf::paraconf
        PDI::pdi
        gslx::paraconfpp
        gslx::geometry_xyvxvy
        gslx::initialisation_xyvxvy
        gslx::interpolation
        gslx::io
        gslx::advection
        gslx::vlasov_xyvxvy
        gslx::poisson_xy
        gslx::quadrature
        gslx::time_integration_xyvxvy
        gslx::utils
)

install(TARGETS landau4d_fft)
