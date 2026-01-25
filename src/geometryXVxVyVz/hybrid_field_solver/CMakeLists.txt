# SPDX-License-Identifier: MIT

add_library("hybrid_field_solver_xvxvyvz" STATIC
    hybridfieldsolver.cpp
    momentscalculator.cpp
    mpimomentscalculator.cpp
)

target_include_directories("hybrid_field_solver_xvxvyvz"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("hybrid_field_solver_xvxvyvz"
    PUBLIC
        DDC::core
        DDC::fft
        gslx::geometry_xvxvyvz
        gslx::pde_solvers
        gslx::quadrature
        gslx::speciesinfo
        gslx::utils
        MPI::MPI_CXX
)

add_library("gslx::hybrid_field_solver_xvxvyvz" ALIAS "hybrid_field_solver_xvxvyvz")
