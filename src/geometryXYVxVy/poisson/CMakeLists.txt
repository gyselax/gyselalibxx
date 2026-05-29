# SPDX-License-Identifier: MIT

foreach(PRECISION_VARIANT IN LISTS PRECISION_VARIANTS_LIST)

add_library("poisson_xy_${PRECISION_VARIANT}" STATIC
    chargedensitycalculator.cpp
    mpichargedensitycalculator.cpp
    nullqnsolver.cpp
    qnsolver.cpp
)

target_include_directories("poisson_xy_${PRECISION_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("poisson_xy_${PRECISION_VARIANT}"
    PUBLIC
        DDC::core
        DDC::fft

        gslx::geometry_xyvxvy_${PRECISION_VARIANT}
        gslx::pde_solvers
        gslx::quadrature
        gslx::speciesinfo
        gslx::utils

        MPI::MPI_CXX
)

add_library("gslx::poisson_xy_${PRECISION_VARIANT}" ALIAS "poisson_xy_${PRECISION_VARIANT}")

endforeach()
