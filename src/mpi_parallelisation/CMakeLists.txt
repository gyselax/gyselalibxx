# SPDX-License-Identifier: MIT

add_library("mpi_parallelisation" INTERFACE)

target_sources("mpi_parallelisation"
    INTERFACE
        FILE_SET HEADERS
            FILES
                impilayout.hpp
                impitranspose.hpp
                mpilayout.hpp
                mpitools.hpp
                mpitransposealltoall.hpp
)

target_link_libraries("mpi_parallelisation"
    INTERFACE
        DDC::core
        gslx::utils
        MPI::MPI_CXX
)

add_library("gslx::mpi_parallelisation" ALIAS "mpi_parallelisation")

install(TARGETS mpi_parallelisation EXPORT GyselalibxxTargets FILE_SET HEADERS)
