# SPDX-License-Identifier: MIT

add_library("initialisation_vparmu" STATIC
    maxwellianequilibrium.cpp
    noperturbinitialization.cpp
)

target_include_directories("initialisation_vparmu"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("initialisation_vparmu"
    PUBLIC
        DDC::core
        gslx::speciesinfo
        gslx::geometry_vparmu
        gslx::utils
)

add_library("gslx::initialisation_vparmu" ALIAS "initialisation_vparmu")
