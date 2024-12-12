# SPDX-License-Identifier: MIT

add_library("advection" INTERFACE)

target_include_directories("advection"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("advection"
    INTERFACE
        DDC::core
        sll::SLL
        gslx::interpolation
        gslx::speciesinfo
        gslx::timestepper
        gslx::utils
)

add_library("gslx::advection" ALIAS "advection")
