# SPDX-License-Identifier: MIT

add_library("initialisation_Kelvin_Helmholtz" INTERFACE)
target_include_directories("initialisation_Kelvin_Helmholtz"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("initialisation_Kelvin_Helmholtz" INTERFACE
    DDC::core
    gslx::geometry_XY
    gslx::utils

)
add_library("gslx::initialisation_Kelvin_Helmholtz" ALIAS "initialisation_Kelvin_Helmholtz")
