# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("fluidrhs_${GEOMETRY_VARIANT}" STATIC
    kinetic_fluid_coupling_source.cpp
)

target_include_directories("fluidrhs_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("fluidrhs_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
        sll::SLL
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::quadrature
        gslx::speciesinfo
        gslx::initialization_${GEOMETRY_VARIANT}
        gslx::reactionrates_${GEOMETRY_VARIANT}
        gslx::timestepper
        gslx::utils_${GEOMETRY_VARIANT}
        gslx::utils
)

add_library("gslx::fluidrhs_${GEOMETRY_VARIANT}" ALIAS "fluidrhs_${GEOMETRY_VARIANT}")

endforeach()
