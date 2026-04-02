# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("poisson_${GEOMETRY_VARIANT}" STATIC
    chargedensitycalculator.cpp
    nullqnsolver.cpp
    qnsolver.cpp
)

target_sources("poisson_${GEOMETRY_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                chargedensitycalculator.hpp
                ichargedensitycalculator.hpp
                iqnsolver.hpp
                nullqnsolver.hpp
                qnsolver.hpp
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

install(TARGETS poisson_${GEOMETRY_VARIANT} EXPORT GyselalibxxTargets FILE_SET HEADERS)

endforeach()
