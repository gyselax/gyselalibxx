# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("initialisation_${GEOMETRY_VARIANT}" STATIC
    maxwellianequilibrium.cpp
    bumpontailequilibrium.cpp
    restartinitialisation.cpp
    singlemodeperturbinitialisation.cpp
)

target_include_directories("initialisation_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
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

endforeach()
