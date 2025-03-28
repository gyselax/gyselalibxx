# SPDX-License-Identifier: MIT

add_executable(landau4d_fft landau4d_fft.cpp)
target_link_libraries(landau4d_fft
    PUBLIC
        DDC::core
        DDC::pdi
        paraconf::paraconf
        PDI::pdi
        gslx::advection
        gslx::geometry_xyvxvy
        gslx::initialisation_xyvxvy
        gslx::interpolation
        gslx::io
        gslx::mpi_parallelisation
        gslx::paraconfpp
        gslx::poisson_xy
        gslx::quadrature
        gslx::time_integration_xyvxvy
        gslx::utils
        gslx::vlasov_xyvxvy
)

install(TARGETS landau4d_fft)
