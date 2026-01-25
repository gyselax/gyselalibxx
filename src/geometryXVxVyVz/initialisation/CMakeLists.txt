# SPDX-License-Identifier: MIT

add_library("initialisation_xvxvyvz" STATIC
    maxwellianequilibrium.cpp
    singlemodeperturbinitialisation.cpp
    hybridmodel_field_initialisation.cpp
)

target_include_directories("initialisation_xvxvyvz"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("initialisation_xvxvyvz"
    PUBLIC
        DDC::core
        gslx::speciesinfo
        gslx::geometry_xvxvyvz
        gslx::utils
)

add_library("gslx::initialisation_xvxvyvz" ALIAS "initialisation_xvxvyvz")
