# SPDX-License-Identifier: MIT

add_library("poisson_xy" STATIC
    chargedensitycalculator.cpp
    mpichargedensitycalculator.cpp
    nullqnsolver.cpp
    qnsolver.cpp
)

target_sources("poisson_xy"
    INTERFACE
        FILE_SET HEADERS
            FILES
                chargedensitycalculator.hpp
                ichargedensitycalculator.hpp
                iqnsolver.hpp
                mpichargedensitycalculator.hpp
                nullqnsolver.hpp
                qnsolver.hpp
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

install(TARGETS poisson_xy EXPORT GyselalibxxTargets FILE_SET HEADERS)
