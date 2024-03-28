# SPDX-License-Identifier: MIT

add_library("collisions_rhs_5D" STATIC
    collisions.cpp
)

target_compile_features("collisions_rhs_5D"
    PUBLIC
        cxx_std_17
)

target_include_directories("collisions_rhs_5D"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("collisions_rhs_5D"
    PUBLIC
        DDC::DDC
        gslx::geometry5D
        gslx::utils
        koliop::koliop
)

add_library("gslx::collisions_rhs_5D" ALIAS "collisions_rhs_5D")
