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

target_sources("rhs_${GEOMETRY_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                collisions_inter.hpp
                collisions_intra.hpp
                collisions_utils.hpp
                irighthandside.hpp
                kinetic_source.hpp
                krook_source_adaptive.hpp
                krook_source_constant.hpp
                mask_tanh.hpp
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

install(TARGETS rhs_${GEOMETRY_VARIANT} EXPORT GyselalibxxTargets FILE_SET HEADERS)

endforeach()
