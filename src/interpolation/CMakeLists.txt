# SPDX-License-Identifier: MIT

add_library("interpolation" INTERFACE)

target_compile_features("interpolation"
    INTERFACE
        cxx_std_17
)

target_include_directories("interpolation"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("interpolation"
    INTERFACE
        DDC::DDC
        sll::splines
        gslx::utils
)

add_library("gslx::interpolation" ALIAS "interpolation")
