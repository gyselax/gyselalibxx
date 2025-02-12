# SPDX-License-Identifier: MIT

add_library("rhs_xperiod_vx" STATIC
    krook_source_adaptive.cpp
    krook_source_constant.cpp
    mask_tanh.cpp
)

target_include_directories("rhs_xperiod_vx"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("rhs_xperiod_vx"
    PUBLIC
        DDC::core
        gslx::geometry_xperiod_vx
        gslx::quadrature
        gslx::speciesinfo
        gslx::initialization_xperiod_vx
        gslx::timestepper
        gslx::utils_xperiod_vx
        gslx::utils
)

add_library("gslx::rhs_xperiod_vx" ALIAS "rhs_xperiod_vx")
