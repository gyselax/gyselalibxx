# SPDX-License-Identifier: MIT
add_library(gyroaverage INTERFACE)

target_include_directories(gyroaverage
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(gyroaverage
    INTERFACE
        DDC::core
        gslx::utils
        gslx::coord_transformations
)

add_library(gslx::gyroaverage ALIAS gyroaverage)

