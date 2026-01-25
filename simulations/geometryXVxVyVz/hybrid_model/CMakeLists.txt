# SPDX-License-Identifier: MIT

add_executable(hybrid_1d3v hybrid_1d3v.cpp)
target_link_libraries(hybrid_1d3v
    PUBLIC
        DDC::core
        DDC::pdi
        paraconf::paraconf
        PDI::pdi
        gslx::advection
        gslx::geometry_xvxvyvz
        gslx::initialisation_xvxvyvz
        gslx::interpolation
        gslx::io
        gslx::mpi_parallelisation
        gslx::paraconfpp
        gslx::hybrid_field_solver_xvxvyvz
        gslx::quadrature
        gslx::time_integration_xvxvyvz
        gslx::utils
        gslx::vlasov_xvxvyvz
)

install(TARGETS hybrid_1d3v)
