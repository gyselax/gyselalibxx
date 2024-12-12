# SPDX-License-Identifier: MIT

add_library("initialization_xyvxvy" STATIC
    maxwellianequilibrium.cpp
    singlemodeperturbinitialization.cpp
)

target_include_directories("initialization_xyvxvy"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("initialization_xyvxvy"
    PUBLIC
        DDC::core
        gslx::speciesinfo
        gslx::geometry_xyvxvy
        gslx::utils
)

add_library("gslx::initialization_xyvxvy" ALIAS "initialization_xyvxvy")
