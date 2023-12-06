# SPDX-License-Identifier: MIT

add_library("advection" INTERFACE)

target_compile_features("advection"
    INTERFACE
        cxx_std_17
)

target_include_directories("advection"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("advection"
    INTERFACE
        DDC::DDC
        sll::splines
        gslx::interpolation
        gslx::speciesinfo
        gslx::utils
)

add_library("gslx::advection" ALIAS "advection")
