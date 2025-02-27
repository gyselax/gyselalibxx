# SPDX-License-Identifier: MIT

add_library("initialisation_xyvxvy" STATIC
    maxwellianequilibrium.cpp
    singlemodeperturbinitialisation.cpp
)

target_include_directories("initialisation_xyvxvy"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("initialisation_xyvxvy"
    PUBLIC
        DDC::core
        gslx::speciesinfo
        gslx::geometry_xyvxvy
        gslx::utils
)

add_library("gslx::initialisation_xyvxvy" ALIAS "initialisation_xyvxvy")
