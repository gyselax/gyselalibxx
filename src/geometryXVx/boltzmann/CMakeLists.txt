# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("boltzmann_${GEOMETRY_VARIANT}" STATIC
    splitvlasovsolver.cpp
    splitrighthandsidesolver.cpp
)

target_compile_features("boltzmann_${GEOMETRY_VARIANT}"
    PUBLIC
        cxx_std_17
)

target_include_directories("boltzmann_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("boltzmann_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        gslx::interpolation
        gslx::speciesinfo
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::rhs_${GEOMETRY_VARIANT}
        gslx::advection
)

add_library("gslx::boltzmann_${GEOMETRY_VARIANT}" ALIAS "boltzmann_${GEOMETRY_VARIANT}")

endforeach()
