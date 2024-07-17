# SPDX-License-Identifier: MIT

add_library("initialization_vparmu" STATIC
    maxwellianequilibrium.cpp
    noperturbinitialization.cpp
)

target_include_directories("initialization_vparmu"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("initialization_vparmu"
    PUBLIC
        DDC::DDC
        gslx::speciesinfo
        gslx::geometry_vparmu
        gslx::utils
)

add_library("gslx::initialization_vparmu" ALIAS "initialization_vparmu")
