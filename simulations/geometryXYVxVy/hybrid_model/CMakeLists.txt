# SPDX-License-Identifier: MIT

add_executable(hybrid4d_fft hybrid4d_fft.cpp)
target_link_libraries(hybrid4d_fft
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
        gslx::hybrid_field_solver_xyvxvy
        gslx::quadrature
        gslx::time_integration_xyvxvy
        gslx::utils
        gslx::vlasov_xyvxvy
)

install(TARGETS hybrid4d_fft)
