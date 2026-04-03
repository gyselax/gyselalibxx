# SPDX-License-Identifier: MIT

add_library("initialisation_Kelvin_Helmholtz" INTERFACE)

target_compile_features("initialisation_Kelvin_Helmholtz" INTERFACE cxx_std_20)

target_sources("initialisation_Kelvin_Helmholtz"
    INTERFACE
        FILE_SET HEADERS
            FILES
                initialisation_Kelvin_Helmholtz.hpp
)

target_link_libraries("initialisation_Kelvin_Helmholtz" INTERFACE
    DDC::core
    gslx::geometry_XY
    gslx::utils

)
add_library("gslx::initialisation_Kelvin_Helmholtz" ALIAS "initialisation_Kelvin_Helmholtz")

install(TARGETS initialisation_Kelvin_Helmholtz EXPORT GyselalibxxTargets FILE_SET HEADERS)
