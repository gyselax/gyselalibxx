# SPDX-License-Identifier: MIT


add_library("interpolation_2D_rp" INTERFACE
)

target_include_directories("interpolation_2D_rp"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("interpolation_2D_rp"
    INTERFACE
        DDC::core
        sll::SLL
        gslx::geometry_RTheta
        gslx::utils

)

add_library("gslx::interpolation_2D_rp" ALIAS "interpolation_2D_rp")
