# SPDX-License-Identifier: MIT


add_library("interpolation_2D_rp" INTERFACE
)

target_compile_features("interpolation_2D_rp"
    INTERFACE
        cxx_std_17
)

target_include_directories("interpolation_2D_rp"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("interpolation_2D_rp"
    INTERFACE
        DDC::DDC
        sll::SLL
        gslx::geometry_RTheta
)

add_library("gslx::interpolation_2D_rp" ALIAS "interpolation_2D_rp")
