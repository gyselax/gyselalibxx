# SPDX-License-Identifier: MIT


add_library("interpolation_RTheta" INTERFACE
)

target_include_directories("interpolation_RTheta"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("interpolation_RTheta"
    INTERFACE
        DDC::core
        gslx::geometry_RTheta
        gslx::utils

)

add_library("gslx::interpolation_RTheta" ALIAS "interpolation_RTheta")
