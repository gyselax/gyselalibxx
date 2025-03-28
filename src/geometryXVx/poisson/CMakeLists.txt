# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("poisson_${GEOMETRY_VARIANT}" STATIC
    chargedensitycalculator.cpp
    nullqnsolver.cpp
    qnsolver.cpp
)

target_include_directories("poisson_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("poisson_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::fft
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::pde_solvers
        gslx::quadrature
        gslx::speciesinfo
        gslx::utils

)

add_library("gslx::poisson_${GEOMETRY_VARIANT}" ALIAS "poisson_${GEOMETRY_VARIANT}")

endforeach()
