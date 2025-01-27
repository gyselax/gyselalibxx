# SPDX-License-Identifier: MIT

add_library("data_types" INTERFACE)

target_include_directories("data_types"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("data_types"
    INTERFACE
        DDC::core
        DDC::splines
        gslx::utils
)

add_library("gslx::data_types" ALIAS "data_types")

