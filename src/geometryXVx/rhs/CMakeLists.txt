# SPDX-License-Identifier: MIT
foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("rhs_${GEOMETRY_VARIANT}" STATIC
    collisions_inter.cpp
    collisions_intra.cpp
    collisions_utils.cpp
    kinetic_source.cpp
    krook_source_adaptive.cpp
    krook_source_constant.cpp
    mask_tanh.cpp
)

target_include_directories("rhs_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("rhs_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::quadrature
        gslx::speciesinfo
        gslx::initialisation_${GEOMETRY_VARIANT}
        gslx::matrix_tools
        gslx::timestepper
        gslx::utils_${GEOMETRY_VARIANT}
        gslx::utils
)

add_library("gslx::rhs_${GEOMETRY_VARIANT}" ALIAS "rhs_${GEOMETRY_VARIANT}")

endforeach()
