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
        sll::splines
        vcx::interpolation_${GEOMETRY_VARIANT}
        vcx::speciesinfo
        vcx::geometry_${GEOMETRY_VARIANT}
        vcx::rhs_${GEOMETRY_VARIANT}
        vcx::advection_${GEOMETRY_VARIANT}
)

add_library("vcx::boltzmann_${GEOMETRY_VARIANT}" ALIAS "boltzmann_${GEOMETRY_VARIANT}")

endforeach()
