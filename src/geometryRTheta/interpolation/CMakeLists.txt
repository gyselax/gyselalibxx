# SPDX-License-Identifier: MIT


add_library("interpolation_rtheta" INTERFACE
)

target_include_directories("interpolation_rtheta"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("interpolation_rtheta"
    INTERFACE
        DDC::core
        gslx::geometry_RTheta
        gslx::utils

)

add_library("gslx::interpolation_rtheta" ALIAS "interpolation_rtheta")
