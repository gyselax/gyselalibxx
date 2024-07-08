# SPDX-License-Identifier: MIT

add_library("initialization_Kelvin_Helmholtz" INTERFACE)
target_include_directories("initialization_Kelvin_Helmholtz"
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries("initialization_Kelvin_Helmholtz" INTERFACE
    DDC::DDC
    gslx::geometry_XY
)
add_library("gslx::initialization_Kelvin_Helmholtz" ALIAS "initialization_Kelvin_Helmholtz")
