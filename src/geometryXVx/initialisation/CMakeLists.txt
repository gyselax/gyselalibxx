# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("initialisation_${GEOMETRY_VARIANT}" STATIC
    iequilibrium.cpp
    maxwellianequilibrium.cpp
    bumpontailequilibrium.cpp
    restartinitialisation.cpp
    singlemodeperturbinitialisation.cpp
)

target_sources("initialisation_${GEOMETRY_VARIANT}"
    INTERFACE
        FILE_SET HEADERS
            FILES
                bumpontailequilibrium.hpp
                iequilibrium.hpp
                iinitialisation.hpp
                maxwellianequilibrium.hpp
                restartinitialisation.hpp
                singlemodeperturbinitialisation.hpp
)

target_link_libraries("initialisation_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::speciesinfo
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::utils
)

add_library("gslx::initialisation_${GEOMETRY_VARIANT}" ALIAS "initialisation_${GEOMETRY_VARIANT}")

install(TARGETS initialisation_${GEOMETRY_VARIANT} EXPORT GyselalibxxTargets FILE_SET HEADERS)

endforeach()
