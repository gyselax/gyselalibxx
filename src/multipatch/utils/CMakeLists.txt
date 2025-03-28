# SPDX-License-Identifier: MIT

add_library("multipatch_utils" INTERFACE)

target_include_directories("multipatch_utils"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("multipatch_utils"
    INTERFACE
        gslx::multipatch_data_types
        gslx::utils
)

add_library("gslx::multipatch_utils" ALIAS "multipatch_utils")


