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

target_compile_features("rhs_${GEOMETRY_VARIANT}"
    PUBLIC
        cxx_std_17
)

target_include_directories("rhs_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("rhs_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        vcx::geometry_${GEOMETRY_VARIANT}
        vcx::quadrature
        vcx::speciesinfo
        vcx::initialization_${GEOMETRY_VARIANT}
        vcx::time_integration_${GEOMETRY_VARIANT}
        vcx::utils_${GEOMETRY_VARIANT}
)

add_library("vcx::rhs_${GEOMETRY_VARIANT}" ALIAS "rhs_${GEOMETRY_VARIANT}")

endforeach()
