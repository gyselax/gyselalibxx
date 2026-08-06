# SPDX-License-Identifier: MIT

add_library("mpi_parallelisation"
  STATIC
    mpi_scope_guard.cpp
)

target_compile_features("mpi_parallelisation" INTERFACE cxx_std_20)

target_sources("mpi_parallelisation"
    INTERFACE
        FILE_SET HEADERS
            FILES
                impilayout.hpp
                impitranspose.hpp
                mpi_scope_guard.hpp
                mpilayout.hpp
                mpitools.hpp
                mpitransposealltoall.hpp
)

target_link_libraries("mpi_parallelisation"
    PUBLIC
        DDC::core
        gyselalibxx::utils
        MPI::MPI_CXX
)

add_library("gyselalibxx::mpi_parallelisation" ALIAS "mpi_parallelisation")

install(TARGETS mpi_parallelisation EXPORT GyselalibxxTargets FILE_SET HEADERS)
