# SPDX-License-Identifier: MIT

add_library("interpolation" INTERFACE)

target_include_directories("interpolation"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("interpolation"
    INTERFACE
        DDC::DDC
        sll::SLL
        gslx::utils
)

add_library("gslx::interpolation" ALIAS "interpolation")
