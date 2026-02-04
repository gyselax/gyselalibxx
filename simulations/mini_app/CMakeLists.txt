# SPDX-License-Identifier: MIT

add_executable(gys_io
    gys_io.cpp
    fluid_moments.cpp
)

target_link_libraries(gys_io
    PUBLIC
        DDC::core
        PDI::pdi
        paraconf::paraconf
        MPI::MPI_CXX
        gslx::io
        gslx::mpi_parallelisation
        gslx::quadrature
        gslx::speciesinfo
        gslx::utils
)

install(TARGETS gys_io)

