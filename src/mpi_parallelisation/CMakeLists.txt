# SPDX-License-Identifier: MIT

add_library("mpi_parallelisation"
  STATIC
    mpi_scope_guard.cpp
)

target_include_directories("mpi_parallelisation"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("mpi_parallelisation"
    PUBLIC
        DDC::core
        gslx::utils
        MPI::MPI_CXX
)

add_library("gslx::mpi_parallelisation" ALIAS "mpi_parallelisation")
