# SPDX-License-Identifier: MIT

add_executable(gys_io
    gys_io.cpp
)

target_link_libraries(gys_io
    PUBLIC
        DDC::core
        PDI::pdi
        paraconf::paraconf
        MPI::MPI_CXX
        gslx::io
        gslx::speciesinfo
        gslx::utils
        gslx::mpi_parallelisation
)

install(TARGETS gys_io)

