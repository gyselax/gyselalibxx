# SPDX-License-Identifier: MIT

add_library("collision" STATIC
    koliop_interface.cpp
)

target_compile_features("collision"
    PUBLIC
        cxx_std_17
)

target_include_directories("collision"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("collision"
    PUBLIC
        DDC::DDC
        gslx::geometry5D
        gslx::utils
        koliop::koliop
)

add_library("gslx::collision" ALIAS "collision")
