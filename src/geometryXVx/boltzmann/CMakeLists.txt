# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("boltzmann_${GEOMETRY_VARIANT}" STATIC
    splitvlasovsolver.cpp
    splitrighthandsidesolver.cpp
)

target_compile_features("boltzmann_${GEOMETRY_VARIANT}" PUBLIC cxx_std_20)

target_sources("boltzmann_${GEOMETRY_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                iboltzmannsolver.hpp
                splitrighthandsidesolver.hpp
                splitvlasovsolver.hpp
)

target_link_libraries("boltzmann_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        gslx::interpolation
        gslx::speciesinfo
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::rhs_${GEOMETRY_VARIANT}
        gslx::advection
)

add_library("gslx::boltzmann_${GEOMETRY_VARIANT}" ALIAS "boltzmann_${GEOMETRY_VARIANT}")

install(TARGETS boltzmann_${GEOMETRY_VARIANT} EXPORT GyselalibxxTargets FILE_SET HEADERS)

endforeach()
