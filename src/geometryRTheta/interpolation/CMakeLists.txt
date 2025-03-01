# SPDX-License-Identifier: MIT


add_library("interpolation_polar" INTERFACE
)

target_include_directories("interpolation_polar"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("interpolation_polar"
    INTERFACE
        DDC::core
        gslx::geometry_RTheta
        gslx::utils

)

add_library("gslx::interpolation_polar" ALIAS "interpolation_polar")
