# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("poisson_${GEOMETRY_VARIANT}" STATIC
    chargedensitycalculator.cpp
    nullqnsolver.cpp
    qnsolver.cpp
)

target_compile_features("poisson_${GEOMETRY_VARIANT}" PUBLIC cxx_std_20)

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
        gyselalibxx::geometry_${GEOMETRY_VARIANT}
        gyselalibxx::pde_solvers
        gyselalibxx::quadrature
        gyselalibxx::speciesinfo
        gyselalibxx::utils

)

add_library("gyselalibxx::poisson_${GEOMETRY_VARIANT}" ALIAS "poisson_${GEOMETRY_VARIANT}")

install(TARGETS poisson_${GEOMETRY_VARIANT} EXPORT gyselalibxxTargets FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/geometryXVx/poisson)

endforeach()
