# SPDX-License-Identifier: MIT

add_library("advection" INTERFACE)

target_include_directories("advection"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("advection"
    INTERFACE
        DDC::core
        gslx::interpolation
        gslx::math_tools
        gslx::speciesinfo
        gslx::timestepper
        gslx::utils
)

add_library("gslx::advection" ALIAS "advection")
