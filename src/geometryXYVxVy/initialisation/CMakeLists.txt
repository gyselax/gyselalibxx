# SPDX-License-Identifier: MIT

add_library("initialisation_xyvxvy" STATIC
    maxwellianequilibrium.cpp
    singlemodeperturbinitialisation.cpp
)

target_sources("initialisation_xyvxvy"
    INTERFACE
        FILE_SET HEADERS
            FILES
                iequilibrium.hpp
                iinitialisation.hpp
                maxwellianequilibrium.hpp
                singlemodeperturbinitialisation.hpp
)

target_link_libraries("initialisation_xyvxvy"
    PUBLIC
        DDC::core
        gslx::speciesinfo
        gslx::geometry_xyvxvy
        gslx::utils
)

add_library("gslx::initialisation_xyvxvy" ALIAS "initialisation_xyvxvy")

install(TARGETS initialisation_xyvxvy EXPORT GyselalibxxTargets FILE_SET HEADERS)
