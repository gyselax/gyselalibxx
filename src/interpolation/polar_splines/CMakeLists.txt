# SPDX-License-Identifier: MIT

add_library("polar_splines" INTERFACE)

target_include_directories("polar_splines"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("polar_splines"
    INTERFACE
        DDC::core
        sll::SLL
        gslx::utils
)

add_library("gslx::polar_splines" ALIAS "polar_splines")
