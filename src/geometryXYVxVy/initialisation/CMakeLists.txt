# SPDX-License-Identifier: MIT

foreach(PRECISION_VARIANT IN LISTS PRECISION_VARIANTS_LIST)

add_library("initialisation_xyvxvy_${PRECISION_VARIANT}" STATIC
    maxwellianequilibrium.cpp
    singlemodeperturbinitialisation.cpp
)

target_include_directories("initialisation_xyvxvy_${PRECISION_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("initialisation_xyvxvy_${PRECISION_VARIANT}"
    PUBLIC
        DDC::core
        gslx::speciesinfo
        gslx::geometry_xyvxvy_${PRECISION_VARIANT}
        gslx::utils
)

add_library("gslx::initialisation_xyvxvy_${PRECISION_VARIANT}" ALIAS "initialisation_xyvxvy_${PRECISION_VARIANT}")

endforeach()
