# SPDX-License-Identifier: MIT

add_library("poisson_xy" STATIC
    chargedensitycalculator.cpp
    mpichargedensitycalculator.cpp
    nullqnsolver.cpp
    qnsolver.cpp
)

target_include_directories("poisson_xy"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("poisson_xy"
    PUBLIC
        DDC::core
        DDC::fft

        gslx::geometry_xyvxvy
        gslx::pde_solvers
        gslx::quadrature
        gslx::speciesinfo
        gslx::utils

        MPI::MPI_CXX
)

add_library("gslx::poisson_xy" ALIAS "poisson_xy")
