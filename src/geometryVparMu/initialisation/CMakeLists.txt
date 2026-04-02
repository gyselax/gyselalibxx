# SPDX-License-Identifier: MIT

add_library("initialisation_vparmu" STATIC
    maxwellianequilibrium.cpp
    noperturbinitialisation.cpp
)

target_sources("initialisation_vparmu"
    INTERFACE
        FILE_SET HEADERS
            FILES
                iequilibrium.hpp
                iinitialisation.hpp
                maxwellianequilibrium.hpp
                noperturbinitialisation.hpp
)

target_link_libraries("initialisation_vparmu"
    PUBLIC
        DDC::core
        gslx::speciesinfo
        gslx::geometry_vparmu
        gslx::utils
)

add_library("gslx::initialisation_vparmu" ALIAS "initialisation_vparmu")

install(TARGETS initialisation_vparmu EXPORT GyselalibxxTargets FILE_SET HEADERS)
