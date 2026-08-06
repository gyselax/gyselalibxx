# SPDX-License-Identifier: MIT

foreach(PRECISION_VARIANT IN LISTS PRECISION_VARIANTS_LIST)

add_library("poisson_xy_${PRECISION_VARIANT}" STATIC
    chargedensitycalculator.cpp
    mpichargedensitycalculator.cpp
    nullqnsolver.cpp
    qnsolver.cpp
)

target_compile_features("poisson_xy_${PRECISION_VARIANT}" PUBLIC cxx_std_20)

target_sources("poisson_xy_${PRECISION_VARIANT}"
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

target_link_libraries("poisson_xy_${PRECISION_VARIANT}"
    PUBLIC
        DDC::core
        DDC::fft

        gyselalibxx::geometry_xyvxvy_${PRECISION_VARIANT}
        gyselalibxx::pde_solvers
        gyselalibxx::quadrature
        gyselalibxx::speciesinfo
        gyselalibxx::utils

        MPI::MPI_CXX
)

add_library("gyselalibxx::poisson_xy_${PRECISION_VARIANT}" ALIAS "poisson_xy_${PRECISION_VARIANT}")

install(TARGETS poisson_xy_${PRECISION_VARIANT} EXPORT GyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryXYVxVy/poisson)

endforeach()
